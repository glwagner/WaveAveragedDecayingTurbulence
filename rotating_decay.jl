using Oceananigans
using Oceananigans.Units
using Statistics
using Printf

arch = GPU()
Nx = Ny = Nz = 256
grid = RectilinearGrid(arch, size=(Nx, Ny, Nz), x=(0, 2π), y=(0, 2π), z=(0, 2π))

model = NonhydrostaticModel(; grid,
                            timestepper = :RungeKutta3,
                            advection = WENO(),
                            coriolis = FPlane(f=1.0))

ϵ(x, y, z) = randn()
set!(model, u=ϵ, v=ϵ, w=ϵ)

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
δ = ∂z(w)

ξ² = Field(ξ^2)
η² = Field(η^2)
ζ² = Field(ζ^2)
δ² = Field(δ^2)
ω² = Field(ξ² + η² + ζ²)
compute!(ω²)

@show ω₀ = sqrt(mean(ω²))
parent(u) ./= ω₀
parent(v) ./= ω₀
parent(w) ./= ω₀

compute!(ω²)
@show sqrt(mean(ω²))

X² = Field(Average(ξ²))
Y² = Field(Average(η²))
Z² = Field(Average(ζ²))
D² = Field(Average(δ²))
U² = Field(Average(u^2))
V² = Field(Average(v^2))
W² = Field(Average(w^2))

u² = Field(u^2)
v² = Field(v^2)
w² = Field(w^2)

e = Field(Average((u² + v² + w²) / 2))

simulation = Simulation(model, Δt=1e-1, stop_time=10000)

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

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

#=
statistics_outputs = (; ξ²=X², η²=Y², ζ²=Z², δ²=D², u²=U², v²=V², w²=W², e)
simulation.output_writers[:statistics] = JLD2OutputWriter(model, statistics_outputs;
                                                          schedule = TimeInterval(0.1),
                                                          filename = "rotating_decay_statistics",
                                                          overwrite_existing = true)
=#

field_outputs = (; u, v, w, ζ, δ)
simulation.output_writers[:slice] = JLD2OutputWriter(model, field_outputs;
                                                     schedule = TimeInterval(10),
                                                     indices = (:, :, Int(Nz/2)),
                                                     filename = "rotating_decay_slice",
                                                     overwrite_existing = true)

simulation.output_writers[:fields] = JLD2OutputWriter(model, field_outputs;
                                                      schedule = TimeInterval(100),
                                                      filename = "rotating_decay_fields",
                                                      overwrite_existing = true)

run!(simulation)

