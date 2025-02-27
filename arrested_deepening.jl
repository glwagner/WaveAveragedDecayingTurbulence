using Oceananigans
using Oceananigans.Units
using Statistics
using Printf
using JLD2

arch = GPU()
Nx = Ny = Nz = 256
x = y = (0, 128)
z = (-64, 0)
topology = (Periodic, Periodic, Bounded)
grid = RectilinearGrid(arch, size=(Nx, Ny, Nz); x, y, z, topology)
weno_order = 9

@inline τ(x, y, t, p) = ifelse(t < p.t★, - p.u★^2, zero(t))
@inline Q(x, y, t, p) = ifelse(t < p.t★, + p.Q₀, zero(t))

# Wave parameters: σ = 2π/8, k ≈ 2π/100, a = 2 meters, ϵ = 0.13.
g = 9.81     # gravitational acceleration
σ = 2π / 8   # frequency
k = σ^2 / g  # wavenumber
c = σ / k    # phase speed

# Medium and strong cases.
# Note that uˢ = ϵ^2 * c
a_med = 1
ϵ_med = a_med * k

a_str = 2
ϵ_str = a_str * k

struct StokesShear
    k :: Float64
    uˢ :: Float64
end

@inline (ss::StokesShear)(z, t) = 2 * ss.k * ss.uˢ * exp(2 * ss.k * z)

cases = [
    #"wind_medium_waves",
    #"wind_strong_waves",
    #"wind_only",
    #"cooling_medium_waves",
    "cooling_strong_waves",
    #"cooling_only",
]

for case in cases

    # Forcing parameters
    N² = 1e-5
    h★ = 24 # some buffer
    w★ = 0.01
    Q★ = w★^3 / h★

    # Note that for wind forcing
    # h★² ~ u★^2 * t★ / N
    # → t★ = h★² * N / u★^2
    #
    # For convection
    # h★² ~ Q * t★ / N²
    # → t★ = h★² * N² / Q

    t★ = 0

    if case == "wind_only"
        stokes_drift = nothing
        u★ = 0.01
        Q₀ = 0.0
        C★ = 1.05
        t★ = h★^2 * sqrt(N²) / (C★ * u★^2)
    elseif case == "cooling_only"
        stokes_drift = nothing
        u★ = 0.0
        Q₀ = Q★
        C★ = 3.0
        t★ = h★^2 * N² / (C★ * Q₀)
    elseif case == "wind_medium_waves"
        stokes_drift = UniformStokesDrift(∂z_uˢ=StokesShear(k, c * ϵ_med^2))
        u★ = 0.01
        Q₀ = 0.0
        C★ = 1.3 # guessing
        t★ = h★^2 * sqrt(N²) / (C★ * u★^2)
    elseif case == "wind_strong_waves"
        stokes_drift = UniformStokesDrift(∂z_uˢ=StokesShear(k, c * ϵ_str^2))
        u★ = 0.01
        Q₀ = 0.0
        C★ = 1.5 # guessing
        t★ = h★^2 * sqrt(N²) / (C★ * u★^2)
    elseif case == "cooling_medium_waves"
        stokes_drift = UniformStokesDrift(∂z_uˢ=StokesShear(k, c * ϵ_med^2))
        u★ = 0.0
        Q₀ = Q★
        C★ = 3.0
        t★ = h★^2 * N² / (C★ * Q₀)
    elseif case == "cooling_strong_waves"
        stokes_drift = UniformStokesDrift(∂z_uˢ=StokesShear(k, c * ϵ_str^2))
        u★ = 0.0
        Q₀ = Q★
        C★ = 3.0
        t★ = h★^2 * N² / (C★ * Q₀)
    end
        
    stop_time = t★ + 48hours
    @info string("Forcing acts until ", prettytime(t★))

    @show parameters = (; t★, u★, Q₀)
    u_top_bc = FluxBoundaryCondition(τ; parameters)
    u_bcs = FieldBoundaryConditions(top=u_top_bc)
    
    b_top_bc = FluxBoundaryCondition(Q; parameters)
    b_bcs = FieldBoundaryConditions(top=b_top_bc)

    # Sponge layer for u, v, w, and T
    gaussian_mask = GaussianMask{:z}(center=-grid.Lz, width=grid.Lz/10)
    u_sponge = v_sponge = w_sponge = Relaxation(rate=1/3minutes, mask=gaussian_mask)

    b_sponge = Relaxation(rate = 1/3minutes,
                          target = LinearTarget{:z}(intercept=0, gradient=N²),
                          mask = gaussian_mask)

    model = NonhydrostaticModel(; grid, stokes_drift,
                                timestepper = :RungeKutta3,
                                tracers = :b,
                                buoyancy = BuoyancyTracer(),
                                boundary_conditions = (; u=u_bcs, b=b_bcs),
                                forcing = (; u=u_sponge, v=v_sponge, w=w_sponge, b=b_sponge),
                                advection = WENO(order=weno_order))

    @show model

    Δz = grid.Lz / grid.Nz
    w★ = max((Δz * Q₀)^(1/3), u★)
    ϵu(x, y, z) = 1e-1 * w★ * rand() * exp(z / 8)
    ϵb(x, y, z) = 1e-1 * N² * Δz * rand() * exp(z / 8)
    bᵢ(x, y, z) = N² * z + ϵb(x, y, z)
    set!(model, u=ϵu, v=ϵu, w=ϵu, b=bᵢ)

    @show model.velocities.u
    @show model.tracers.b

    simulation = Simulation(model, Δt=1.0; stop_time)
    wizard = TimeStepWizard(cfl=0.5, max_change=1.1)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    start_time = Ref(time_ns())
    function progress(sim)
        msg = @sprintf("Iter: %d, time: %s", iteration(sim), prettytime(sim))
        msg *= @sprintf(", Δt: %s", prettytime(sim.Δt))

        u, v, w = sim.model.velocities
        msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e)",
                        maximum(abs, interior(u)),
                        maximum(abs, interior(v)),
                        maximum(abs, interior(w)))

        elapsed = (time_ns() - start_time[]) / 1e9
        msg *= @sprintf(", wall time: %s", prettytime(elapsed))
        start_time[] = time_ns()

        @info msg

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    prefix = "arrested_deepening_$(case)"

    u, v, w = model.velocities

    η = ∂x(w) - ∂z(u)
    b = model.tracers.b
    outputs = (; u, v, w, b, η)

    U = Average(u, dims=(1, 2))
    B = Average(b, dims=(1, 2))
    W² = Average(w^2, dims=(1, 2))
    u′ = Field(u - Field(U))
    e = Field((u′^2 + v^2 + w^2) / 2)
    E = Average(e, dims=(1, 2))

    slices_time_interval = 10
    fields_time_interval = 100

    function init(file, model)
        file["u★"] = u★
        file["t★"] = t★
        file["Q₀"] = Q₀
    end

    save_interval = 2minutes

    simulation.output_writers[:averages] = JLD2OutputWriter(model, (u=U, b=B, w²=W², e=E); init,
                                                           schedule = TimeInterval(save_interval),
                                                           filename = prefix * "_statistics",
                                                           overwrite_existing = true)

    simulation.output_writers[:xy] = JLD2OutputWriter(model, outputs; init,
                                                      schedule = TimeInterval(save_interval),
                                                      indices = (:, :, Nz),
                                                      filename = prefix * "_xy",
                                                      overwrite_existing = true)

    simulation.output_writers[:xz] = JLD2OutputWriter(model, outputs; init,
                                                      schedule = TimeInterval(save_interval),
                                                      indices = (:, 1, :),
                                                      filename = prefix * "_xz",
                                                      overwrite_existing = true)

    simulation.output_writers[:yz] = JLD2OutputWriter(model, outputs; init,
                                                      schedule = TimeInterval(save_interval),
                                                      indices = (1, :, :),
                                                      filename = prefix * "_yz",
                                                      overwrite_existing = true)

    @info "Running $case..."

    run!(simulation)
end

