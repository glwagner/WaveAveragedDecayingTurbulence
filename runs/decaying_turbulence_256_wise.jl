using FourierFlows
using FFTW
using Oceananigans
using Oceananigans.Units
using Statistics
using Printf
using JLD2

using Oceananigans.Operators

@inline ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2

@inline function turbulent_kinetic_energyᶜᶜᶜ(i, j, k, grid, u, v, w)
    u² = ℑxᶜᵃᵃ(i, j, k, grid, ϕ², u)
    v² = ℑyᵃᶜᵃ(i, j, k, grid, ϕ², v)
    w² = ℑzᵃᵃᶜ(i, j, k, grid, ϕ², w)
    return (u² + v² + w²) / 2
end

@inline x_vorticityᶜᶠᶠ(i, j, k, grid, u, v, w) = ∂zᶜᶠᶠ(i, j, k, grid, v) - ∂yᶜᶠᶠ(i, j, k, grid, w)
@inline y_vorticityᶠᶜᶠ(i, j, k, grid, u, v, w) = ∂xᶠᶜᶠ(i, j, k, grid, w) - ∂zᶠᶜᶠ(i, j, k, grid, u)
@inline z_vorticityᶠᶠᶜ(i, j, k, grid, u, v, w) = ∂xᶠᶠᶜ(i, j, k, grid, v) - ∂yᶠᶠᶜ(i, j, k, grid, u)

@inline x_enstrophyᶜᶠᶠ(i, j, k, grid, u, v, w) = (∂zᶜᶠᶠ(i, j, k, grid, v) - ∂yᶜᶠᶠ(i, j, k, grid, w))^2
@inline y_enstrophyᶠᶜᶠ(i, j, k, grid, u, v, w) = (∂xᶠᶜᶠ(i, j, k, grid, w) - ∂zᶠᶜᶠ(i, j, k, grid, u))^2
@inline z_enstrophyᶠᶠᶜ(i, j, k, grid, u, v, w) = (∂xᶠᶠᶜ(i, j, k, grid, v) - ∂yᶠᶠᶜ(i, j, k, grid, u))^2

@inline function ω²ᶠᶠᶠ(i, j, k, grid, u, v, w)
    ξ² = ℑxᶠᵃᵃ(i, j, k, grid, x_enstrophyᶜᶠᶠ, u, v, w)
    η² = ℑyᵃᶠᵃ(i, j, k, grid, y_enstrophyᶠᶜᶠ, u, v, w)
    ζ² = ℑzᵃᵃᶠ(i, j, k, grid, z_enstrophyᶠᶠᶜ, u, v, w)
    return ξ² + η² + ζ²
end

function set_zero_mean_velocity_and_rms_vorticity!(model, ω_rms=1)
    u, v, w = model.velocities
    U = mean(u); V = mean(v); W = mean(w)
    parent(u) .-= U; parent(v) .-= V; parent(w) .-= W

    ω² = KernelFunctionOperation{Face, Face, Face}(ω²ᶠᶠᶠ, grid, u, v, w)
    ω² = compute!(Field(ω²))
    ω₀ = sqrt(mean(ω²))
    parent(u) .*= ω_rms / ω₀
    parent(v) .*= ω_rms / ω₀
    parent(w) .*= ω_rms / ω₀

    compute!(ω²)
    @show sqrt(mean(ω²))
    return nothing
end

function velocity_spectral_shape(k)
    k₀ = 32 * 2π
    k′ = k / k₀
    return k * sqrt(exp(-2 * k′^2))
end

#####
##### WISE-conference cut-down: only the two cases we need for the (kx,kz) spectrum plot.
#####
kinds = [
    "medium_surface_waves",
]

# Allow overriding stop_time / grid size via env vars for sanity checks.
const STOP_TIME = parse(Float64, get(ENV, "STOP_TIME", "1e4"))
const N         = parse(Int,     get(ENV, "N", "256"))
const SANITY    = get(ENV, "SANITY", "0") == "1"

Nx = Ny = Nz = N
weno_order = 9
cfl = 0.5
arch = Oceananigans.GPU()
stop_time = STOP_TIME

statistics_saves = 1000
n_statistics = range(-2, stop=4, length=statistics_saves)
statistics_times = 10 .^ n_statistics

slices_saves = 400
n_slices = range(-1, stop=4, length=slices_saves)
slices_times = 10 .^ n_slices

initialization_stop_time = 1e3
initialization_stop_ω_rms = 10

timestepper = :RungeKutta3
advection = WENO(order=weno_order)
x = y = z = (0, 1)
topology = (Periodic, Periodic, Bounded)
grid = RectilinearGrid(arch, size=(Nx, Ny, Nz), halo=(7, 7, 7); x, y, z, topology)

initial_conditions_filename = "initial_conditions_$Nx.jld2"

if !isfile(initial_conditions_filename)
    spectral_grid = ThreeDGrid(nx=Nx, Lx=1)
    k = sqrt.(spectral_grid.Krsq)

    FT = eltype(grid)
    θu = randn(Complex{FT}, size(spectral_grid.Krsq))
    θv = randn(Complex{FT}, size(spectral_grid.Krsq))
    θw = randn(Complex{FT}, size(spectral_grid.Krsq))

    û = @. θu * velocity_spectral_shape(k)
    v̂ = @. θv * velocity_spectral_shape(k)
    ŵ = @. θw * velocity_spectral_shape(k)
    û[1, 1, 1] = 0; v̂[1, 1, 1] = 0; ŵ[1, 1, 1] = 0

    u₀  = irfft(û, spectral_grid.nx)
    v₀  = irfft(v̂, spectral_grid.nx)
    w₀′ = irfft(ŵ, spectral_grid.nx)

    w₀ = zeros(Nx, Ny, Nz+1)
    w₀[:, :, 2:Nz] .= w₀′[:, :, 2:Nz]

    init_model = NonhydrostaticModel(; grid, timestepper, advection)
    set!(init_model, u=u₀, v=v₀, w=w₀)
    set_zero_mean_velocity_and_rms_vorticity!(init_model, 1000)

    max_u = maximum(abs, init_model.velocities.u)
    Δx = 1 / Nx
    Δt = 1e-3 * Δx / max_u

    simulation = Simulation(init_model; Δt, stop_time=initialization_stop_time)
    conjure_time_step_wizard!(simulation, IterationInterval(3); cfl)

    start_time = Ref(time_ns())
    u, v, w = init_model.velocities
    ω² = KernelFunctionOperation{Face, Face, Face}(ω²ᶠᶠᶠ, grid, u, v, w)

    function progress(sim)
        msg = @sprintf("Iter: % 6d, time: %7.2f, Δt: %.2e", iteration(sim), time(sim), sim.Δt)
        u, v, w = sim.model.velocities
        msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e)", maximum(abs, u), maximum(abs, v), maximum(abs, w))
        Ω = sqrt(mean(ω²))
        msg *= @sprintf(", rms(ω): %.2e", Ω)
        elapsed = (time_ns() - start_time[]) / 1e9
        msg *= @sprintf(", wall time: %s", prettytime(elapsed))
        start_time[] = time_ns()
        @info msg
        return nothing
    end

    add_callback!(simulation, progress, IterationInterval(100))

    function stop_simulation(sim)
        Ω = sqrt(mean(ω²))
        if Ω < initialization_stop_ω_rms
            @info "rms(ω) < initialization_stop_ω_rms, stopping simulation."
            progress(sim)
            simulation.running = false
        end
        return nothing
    end

    add_callback!(simulation, stop_simulation, IterationInterval(10))

    @info "Generating initial conditions for Nx = $(grid.Nx)..."
    run!(simulation)

    u, v, w = init_model.velocities
    ui = Array(interior(u))
    vi = Array(interior(v))
    wi = Array(interior(w))

    file = jldopen(initial_conditions_filename, "a+")
    file["ui"] = ui; file["vi"] = vi; file["wi"] = wi
    close(file)
else
    file = jldopen(initial_conditions_filename)
    ui = file["ui"]; vi = file["vi"]; wi = file["wi"]
    close(file)
end

struct ShallowStokesShear
    shear :: Float64
end

@inline (s::ShallowStokesShear)(z, t) = s.shear * z

for kind in kinds

    local max_u, ω², Δx, Δt, u, v, w, start_time, simulation

    if kind == "isotropic"
        kwargs = NamedTuple()
    elseif kind == "medium_surface_waves"
        kwargs = (; stokes_drift = UniformStokesDrift(∂z_uˢ=ShallowStokesShear(0.5)))
    end

    model = NonhydrostaticModel(; grid, timestepper, advection, kwargs...)
    set!(model, u=ui, v=vi, w=wi)

    max_u = maximum(abs, model.velocities.u)
    Δx = 1 / Nx
    Δt = 1e-3 * Δx / max_u
    simulation = Simulation(model; Δt, stop_time)
    conjure_time_step_wizard!(simulation, IterationInterval(3); cfl)

    prefix = "decaying_turbulence_$(Nx)_$(weno_order)_$(kind)"

    u, v, w = model.velocities
    η = ∂x(w) - ∂z(u)
    η² = KernelFunctionOperation{Face, Center, Face}(y_enstrophyᶠᶜᶠ, grid, u, v, w)
    ω² = KernelFunctionOperation{Face, Face, Face}(ω²ᶠᶠᶠ, grid, u, v, w)
    e  = KernelFunctionOperation{Center, Center, Center}(turbulent_kinetic_energyᶜᶜᶜ, grid, u, v, w)

    e  = Field(e); η² = Field(η²); ω² = Field(ω²)

    max_η²(model) = maximum(η²)
    max_e(model)  = maximum(e)

    Y² = Field(Average(η²)); Ω² = Field(Average(ω²)); E = Field(Average(e))

    start_time = Ref(time_ns())

    function progress(sim)
        msg = @sprintf("Iter: % 6d, time: %7.2f, Δt: %.2e", iteration(sim), time(sim), sim.Δt)
        u, v, w = sim.model.velocities
        msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e)", maximum(abs, u), maximum(abs, v), maximum(abs, w))
        compute!(e); compute!(ω²)
        E = mean(e); Ω = sqrt(mean(ω²))
        msg *= @sprintf(", ⟨e⟩: %.2e, rms(ω): %.2e", E, Ω)
        elapsed = (time_ns() - start_time[]) / 1e9
        msg *= @sprintf(", wall time: %s", prettytime(elapsed))
        start_time[] = time_ns()
        @info msg
        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

    statistics_outputs = (; e=E, η²=Y², ω²=Ω², max_η², max_e)
    slice_outputs = (; η, u, v, w, e)

    simulation.output_writers[:statistics] = JLD2OutputWriter(model, statistics_outputs;
                                                              schedule = SpecifiedTimes(statistics_times),
                                                              filename = prefix * "_statistics",
                                                              overwrite_existing = true)

    simulation.output_writers[:xz] = JLD2OutputWriter(model, slice_outputs;
                                                     indices = (:, 1, :),
                                                     schedule = SpecifiedTimes(slices_times),
                                                     filename = prefix * "_xz",
                                                     overwrite_existing = true)

    simulation.output_writers[:yz] = JLD2OutputWriter(model, slice_outputs;
                                                     indices = (1, :, :),
                                                     schedule = SpecifiedTimes(slices_times),
                                                     filename = prefix * "_yz",
                                                     overwrite_existing = true)

    simulation.output_writers[:xy] = JLD2OutputWriter(model, slice_outputs;
                                                     indices = (:, :, Nz),
                                                     schedule = SpecifiedTimes(slices_times),
                                                     filename = prefix * "_xy",
                                                     overwrite_existing = true)

    U_avg = Average(u, dims=(1, 2))
    V_avg = Average(v, dims=(1, 2))
    simulation.output_writers[:avg] = JLD2OutputWriter(model, (u=U_avg, v=V_avg);
                                                      schedule = SpecifiedTimes(slices_times),
                                                      filename = prefix * "_averages",
                                                      overwrite_existing = true)

    # Full 3D fields at the times we need for (kx,ky,kz) spectra.
    field_outputs = (; u, v, w)
    field_times = SANITY ? SpecifiedTimes(0.5) : SpecifiedTimes(100, 1000, 10000)
    simulation.output_writers[:fields] = JLD2OutputWriter(model, field_outputs;
                                                          schedule = field_times,
                                                          filename = prefix * "_fields",
                                                          with_halos = true,
                                                          overwrite_existing = true)

    @info "Running $kind turbulence with N = $(grid.Nx), stop_time = $stop_time..."
    run!(simulation)
end
