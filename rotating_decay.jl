using Oceananigans
using Oceananigans.Units
using Statistics
using Printf
using GLMakie

arch = CPU()
Nx = Ny = Nz = 16
grid = RectilinearGrid(arch, size=(Nx, Ny, Nz), x=(0, 2π), y=(0, 2π), z=(0, 2π))

model = NonhydrostaticModel(; grid,
                            timestepper = :RungeKutta3,
                            advection = WENO(),
                            coriolis = FPlane(f=1e-1))

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

ω₀ = sqrt(mean(ξ^2 + η^2 + ζ^2))
parent(u) ./= ω₀
parent(v) ./= ω₀
parent(w) ./= ω₀

@show ω = sqrt(mean(ξ^2 + η^2 + ζ^2))

ξ² = Field(Average(ξ^2))
η² = Field(Average(η^2))
ζ² = Field(Average(ζ^2))
u² = Field(Average(u^2))
v² = Field(Average(v^2))
w² = Field(Average(w^2))

e = Field(Average((u^2 + v^2 + w^2) / 2))

simulation = Simulation(model, Δt=1e-2, stop_time=100)

wizard = TimeStepWizard(cfl=0.7, max_change=1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

function progress(sim)
    msg = @sprintf("Iter: %d, time: %.2f", iteration(sim), time(sim))
    msg *= @sprintf(", Δt: %.4f", sim.Δt)
    @info msg
    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

statistics_output = (; ξ², η², ζ², u², v², w², e)
simulation.output_writers[:statistics] = JLD2OutputWriter(model, statistics_output;
                                                          schedule = TimeInterval(0.1),
                                                          filename = "rotating_decay_statistics",
                                                          overwrite_existing = true)

simulation.output_writers[:slice] = JLD2OutputWriter(model, model.velocities;
                                                     schedule = TimeInterval(0.1),
                                                     indices = (:, :, Nz),
                                                     filename = "rotating_decay_slice",
                                                     overwrite_existing = true)

field_output = (; u, v, w, ζ)
simulation.output_writers[:fields] = JLD2OutputWriter(model, model.velocities;
                                                      schedule = TimeInterval(1),
                                                      filename = "rotating_decay_fields",
                                                      overwrite_existing = true)

run!(simulation)

fig = Figure()
ut = FieldTimeSeries("rotating_decay_fields.jld2", "u")
t = ut.times
Nt = length(t)

axu = Axis(fig[1, 1])
slider = Slider(fig[2, 1], range=1:Nt, startvalue=1)
n = slider.value

un = @lift interior(ut[$n], :, :, 1)

heatmap!(axu, un)

display(fig)

fig = Figure()
axω = Axis(fig[1, 1], yscale=log10)
axe = Axis(fig[2, 1], yscale=log10)

ξ²t = FieldTimeSeries("rotating_decay_statistics.jld2", "ξ²")
η²t = FieldTimeSeries("rotating_decay_statistics.jld2", "η²")
ζ²t = FieldTimeSeries("rotating_decay_statistics.jld2", "ζ²")

u²t = FieldTimeSeries("rotating_decay_statistics.jld2", "u²")
v²t = FieldTimeSeries("rotating_decay_statistics.jld2", "v²")
w²t = FieldTimeSeries("rotating_decay_statistics.jld2", "w²")
et = FieldTimeSeries("rotating_decay_statistics.jld2", "e")

t = et.times

lines!(axω, t, sqrt.(ξ²t[1, 1, 1, :]), linewidth=4, label="x")
lines!(axω, t, sqrt.(η²t[1, 1, 1, :]), linewidth=4, label="y")
lines!(axω, t, sqrt.(ζ²t[1, 1, 1, :]), linewidth=4, label="z")
axislegend(axω)

lines!(axe, t, u²t[1, 1, 1, :] .* 3/2, label="3u²/2", linewidth=4)
lines!(axe, t, v²t[1, 1, 1, :] .* 3/2, label="3v²/2", linewidth=4)
lines!(axe, t, w²t[1, 1, 1, :] .* 3/2, label="3w²/2", linewidth=4)
lines!(axe, t, et[1, 1, 1, :], label="e", linewidth=4)
axislegend(axe)

display(fig)
