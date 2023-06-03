using Oceananigans
using Oceananigans.Units
using Statistics
using Printf
using JLD2

arch = GPU()
Nx = Ny = Nz = 256
x = y = (0, 48)
z = (-24, 0)
topology = (Periodic, Periodic, Bounded)
grid = RectilinearGrid(arch, size=(Nx, Ny, Nz); x, y, z, topology)

N² = 1e-4
u★ = 0.01
u_top_bc = FluxBoundaryCondition(-u★^2)
u_bcs = FieldBoundaryConditions(top=u_top_bc)

h★ = 12
# h★² ~ u★^2 * t★ / N
# → t★ = h★^2 * N / u★^2
t★ = h★^2 * sqrt(N²) / u★^2
@show prettytime(t★)

uˢ = 1.0
g = 9.81
σm = 2π / 12
σl = 2π / 24

k_medium = σm^2 / g
k_long   = σl^2 / g

@show h_medium = 1 / 2k_medium
@show h_long   = 1 / 2k_long

struct StokesShear
    k :: Float64
    uˢ :: Float64
end

@inline (ss::StokesShear)(z, t) = 2 * ss.k * ss.uˢ * exp(2 * ss.k * z)

uˢ_strong = uˢ * k_medium / k_long

@show uˢ_strong
@show uˢ_weak

cases = [
    "medium",
    "long",
    "strong",
    "none",
]

for case in cases

    if case == "medium"
        stokes_drift = UniformStokesDrift(∂z_uˢ=StokesShear(k_medium, uˢ))
    elseif case == "long"
        stokes_drift = UniformStokesDrift(∂z_uˢ=StokesShear(k_long, uˢ))
    elseif case == "strong"
        stokes_drift = UniformStokesDrift(∂z_uˢ=StokesShear(k_long, uˢ_strong))
    elseif case == "none"
        stokes_drift = nothing
    end

    model = NonhydrostaticModel(; grid, stokes_drift,
                                timestepper = :RungeKutta3,
                                tracers = :b,
                                buoyancy = BuoyancyTracer(),
                                boundary_conditions = (; u=u_bcs),
                                advection = WENO())

    ϵ(x, y, z) = 1e-2 * u★ * rand() * exp(z / 4)
    bᵢ(x, y, z) = N² * z
    set!(model, u=ϵ, v=ϵ, w=ϵ, b=bᵢ)

    u, v, w = model.velocities
    η = ∂x(w) - ∂z(u)
    e = (u^2 + v^2 + w^2) / 2

    simulation = Simulation(model, Δt=1.0, stop_time=t★)
    wizard = TimeStepWizard(cfl=0.5, max_change=1.1)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    start_time = Ref(time_ns())
    function progress(sim)
        msg = @sprintf("Iter: %d, time: %s", iteration(sim), prettytime(sim))
        msg *= @sprintf(", Δt: %s", prettytime(sim.Δt))

        u, v, w = sim.model.velocities
        msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e)", maximum(abs, u), maximum(abs, v), maximum(abs, w))

        elapsed = (time_ns() - start_time[]) / 1e9
        msg *= @sprintf(", wall time: %s", prettytime(elapsed))
        start_time[] = time_ns()

        @info msg

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    prefix = "wavy_kato_phillips_$(case)_waves"

    b = model.tracers.b
    xy_outputs = (; u, v, w)
    xz_outputs = (; u, v, w)
    field_outputs = (; u, v, w)

    U = Average(u, dims=(1, 2))
    B = Average(b, dims=(1, 2))

    slices_time_interval = 10
    fields_time_interval = 100

    simulation.output_writers[:averages] = JLD2OutputWriter(model, (u=U, b=B);
                                                           schedule = TimeInterval(10minutes),
                                                           filename = prefix * "_statistics",
                                                           overwrite_existing = true)

    simulation.output_writers[:xy] = JLD2OutputWriter(model, merge(model.velocities, model.tracers);
                                                      schedule = TimeInterval(10minutes),
                                                      indices = (:, :, Nz),
                                                      filename = prefix * "_xy",
                                                      overwrite_existing = true)

    simulation.output_writers[:xz] = JLD2OutputWriter(model, merge(model.velocities, model.tracers);
                                                      schedule = TimeInterval(10minutes),
                                                      indices = (:, 1, :),
                                                      filename = prefix * "_xz",
                                                      overwrite_existing = true)

    simulation.output_writers[:yz] = JLD2OutputWriter(model, merge(model.velocities, model.tracers);
                                                      schedule = TimeInterval(10minutes),
                                                      indices = (1, :, :),
                                                      filename = prefix * "_yz",
                                                      overwrite_existing = true)

    @info "Running $case..."

    run!(simulation)
end

