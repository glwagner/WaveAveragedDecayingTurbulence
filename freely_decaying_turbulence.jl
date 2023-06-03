using FourierFlows
using FFTW
using Oceananigans
using Oceananigans.Units
using Statistics
using Printf
using JLD2

arch = Oceananigans.GPU()
Nx = Ny = Nz = 256
x = y = (0, 1)
z = (-1, 0)
topology = (Periodic, Periodic, Bounded)
grid = RectilinearGrid(arch, size=(Nx, Ny, Nz); x, y, z, topology)

# Initial conditions
initial_conditions_filename = "initial_conditions_$Nx.jld2"

rm(initial_conditions_filename, force=true)

if !isfile(initial_conditions_filename)
    Nx, Ny, Nz = size(grid)

    spectral_grid = ThreeDGrid(nx=Nx, Lx=1)
    k = sqrt.(spectral_grid.Krsq)

    FT = eltype(grid)
    θu = randn(Complex{FT}, size(spectral_grid.Krsq))
    θv = randn(Complex{FT}, size(spectral_grid.Krsq))
    θw = randn(Complex{FT}, size(spectral_grid.Krsq))

    û = θu ./ k
    v̂ = θv ./ k
    ŵ = θw ./ k

    @show size(k)

    û[1, 1, 1] = 0
    v̂[1, 1, 1] = 0
    ŵ[1, 1, 1] = 0

    ui  = irfft(û, spectral_grid.nx)
    vi  = irfft(v̂, spectral_grid.nx)
    wi′ = irfft(ŵ, spectral_grid.nx)

    wi = zeros(Nx, Ny, Nz+1)
    wi[:, :, 2:Nz] .= wi′[:, :, 2:Nz]

    @show size(ui)
    @show size(vi)
    @show size(wi)
    @show mean(ui)
    @show mean(vi)
    @show mean(wi)

    file = jldopen(initial_conditions_filename, "a+")
    file["ui"] = ui
    file["vi"] = vi
    file["wi"] = wi
    close(file)

else
    file = jldopen(initial_conditions_filename)
    ui = file["ui"]
    vi = file["vi"]
    wi = file["wi"]
    close(file)
end

function set_zero_mean_velocity_and_unity_vorticity!(model)
    # Subtract mean
    u, v, w = model.velocities
    U = mean(u)
    V = mean(v)
    W = mean(w)
    parent(u) .-= U
    parent(v) .-= V
    parent(w) .-= W

    # Scale rms vorticity
    ξ = ∂z(v) - ∂y(w)
    η = ∂x(w) - ∂z(u)
    ζ = ∂x(v) - ∂y(u)

    ξ² = Field(ξ^2)
    η² = Field(η^2)
    ζ² = Field(ζ^2)
    ω² = Field(ξ² + η² + ζ²)
    compute!(ω²)

    ω₀ = sqrt(mean(ω²))
    parent(u) ./= ω₀
    parent(v) ./= ω₀
    parent(w) ./= ω₀

    compute!(ω²)
    @show sqrt(mean(ω²))

    return nothing
end

coriolis = FPlane(f=0.1)

@inline navid_s(z, t) = 0.2 * (z + 1/2)
navid_stokes_drift = UniformStokesDrift(∂z_uˢ=navid_s)

@inline s(z, t) = 0.1 * (z + 1)
stokes_drift = UniformStokesDrift(∂z_uˢ=s)

@inline weak_s(z, t) = 0.05 * (z + 1)
weak_stokes_drift = UniformStokesDrift(∂z_uˢ=weak_s)

@inline very_weak_s(z, t) = 0.02 * (z + 1)
very_weak_stokes_drift = UniformStokesDrift(∂z_uˢ=very_weak_s)

@inline strong_s(z, t) = 0.2 * (z + 1)
strong_stokes_drift = UniformStokesDrift(∂z_uˢ=strong_s)

@inline very_strong_s(z, t) = 0.5 * (z + 1)
very_strong_stokes_drift = UniformStokesDrift(∂z_uˢ=very_strong_s)

@inline very_deep_s(z, t) = 0.8 * exp(8z)
very_deep_stokes_drift = UniformStokesDrift(∂z_uˢ=very_deep_s)

@inline deep_s(z, t) = 0.4 * exp(4z)
deep_stokes_drift = UniformStokesDrift(∂z_uˢ=deep_s)

kinds = ["isotropic",
         "rotating",
         "surface_waves",
         #"navid_surface_waves",
         "weak_surface_waves",
         "deep_surface_waves",
         "very_deep_surface_waves",
         "strong_surface_waves",
         "very_strong_surface_waves",
]

for kind in kinds

    if kind == "rotating"
        kwargs = (; coriolis)
    elseif kind == "isotropic"
        kwargs = NamedTuple()
    elseif kind == "surface_waves"
        kwargs = (; stokes_drift)
    elseif kind == "navid_surface_waves"
        kwargs = (; stokes_drift=navid_stokes_drift)
    elseif kind == "strong_surface_waves"
        kwargs = (; stokes_drift=strong_stokes_drift)
    elseif kind == "very_strong_surface_waves"
        kwargs = (; stokes_drift=very_strong_stokes_drift)
    elseif kind == "weak_surface_waves"
        kwargs = (; stokes_drift=weak_stokes_drift)
    elseif kind == "very_weak_surface_waves"
        kwargs = (; stokes_drift=very_weak_stokes_drift)
    elseif kind == "deep_surface_waves"
        kwargs = (; stokes_drift=deep_stokes_drift)
    elseif kind == "very_deep_surface_waves"
        kwargs = (; stokes_drift=very_deep_stokes_drift)
    end
                                
    model = NonhydrostaticModel(; grid,
                                timestepper = :RungeKutta3,
                                advection = WENO(),
                                kwargs...)

    set!(model, u=ui, v=vi, w=wi)
    set_zero_mean_velocity_and_unity_vorticity!(model)
    u, v, w = model.velocities
    ξ = ∂z(v) - ∂y(w)
    η = ∂x(w) - ∂z(u)
    ζ = ∂x(v) - ∂y(u)
    δ = ∂z(w)
    e = (u^2 + v^2 + w^2) / 2
    α = ∂y(v)
    δ² = Field(δ^2)
    α² = Field(α^2)
    E = Field(Average(Field(e)))

    simulation = Simulation(model, Δt=1e-2, stop_time=2e4)
    wizard = TimeStepWizard(cfl=0.7, max_change=1.1)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(3))

    start_time = Ref(time_ns())
    function progress(sim)
        msg = @sprintf("Iter: %d, time: %.2f", iteration(sim), time(sim))
        msg *= @sprintf(", Δt: %.4f", sim.Δt)

        u, v, w = sim.model.velocities
        msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e)", maximum(abs, u), maximum(abs, v), maximum(abs, w))

        u, v, w = sim.model.velocities
        msg *= @sprintf(", max|ω|: (%.2e, %.2e, %.2e)", maximum(abs, ξ), maximum(abs, η), maximum(abs, ζ))

        elapsed = (time_ns() - start_time[]) / 1e9
        msg *= @sprintf(", wall time: %s", prettytime(elapsed))
        start_time[] = time_ns()

        @info msg

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

    prefix = "decaying_turbulence_$(Nx)_$(kind)"

    η² = Field(η^2)
    ζ² = Field(ζ^2)

    Y² = Field(Average(η²))
    Z² = Field(Average(ζ²))

    if kind == "rotating"
        max_ζ(model) = maximum(abs, ζ)
        statistics_outputs = (; e=E, ζ²=Z², max_ζ)
        slice_outputs = (; ζ, u, v, w, e)
        indices = (:, :, Int(Nz/2))
        section_indices = (:, 1, :)
    else
        max_η(model) = maximum(abs, η)
        statistics_outputs = (; e=E, η²=Y², max_η)
        slice_outputs = (; η, u, v, w, e)
        indices = (:, 1, :)
        section_indices = (1, :, :)
    end

    statistics_time_interval = 10
    slices_time_interval = 10

    simulation.output_writers[:statistics] = JLD2OutputWriter(model, statistics_outputs;
                                                              schedule = TimeInterval(statistics_time_interval),
                                                              filename = prefix * "_statistics",
                                                              overwrite_existing = true)

    simulation.output_writers[:xy] = JLD2OutputWriter(model, slice_outputs;
                                                      indices = (:, :, Nz),
                                                      schedule = TimeInterval(slices_time_interval),
                                                      filename = prefix * "_xy",
                                                      overwrite_existing = true)

    simulation.output_writers[:xz] = JLD2OutputWriter(model, slice_outputs;
                                                      indices = (:, 1, :),
                                                      schedule = TimeInterval(slices_time_interval),
                                                      filename = prefix * "_xz",
                                                      overwrite_existing = true)

    simulation.output_writers[:yz] = JLD2OutputWriter(model, slice_outputs;
                                                      indices = (1, :, :),
                                                      schedule = TimeInterval(slices_time_interval),
                                                      filename = prefix * "_yz",
                                                      overwrite_existing = true)

    U = Average(u, dims=(1, 2))
    V = Average(v, dims=(1, 2))
    simulation.output_writers[:avg] = JLD2OutputWriter(model, (u=U, v=V);
                                                       schedule = TimeInterval(slices_time_interval),
                                                       filename = prefix * "_averages",
                                                       overwrite_existing = true)

    #=
    field_outputs = (; u, v, w, ζ, η)
    fields_time_interval = 100
    simulation.output_writers[:fields] = JLD2OutputWriter(model, field_outputs;
                                                          schedule = TimeInterval(fields_time_interval),
                                                          filename = prefix * "_fields",
                                                          overwrite_existing = true)
    =#

    @info "Running $kind turbulence..."
    run!(simulation)
end

