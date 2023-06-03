using GLMakie
using Oceananigans
using DifferentialEquations

filename_waves        = "decaying_turbulence_256_surface_waves_statistics.jld2"
filename_strong_waves = "decaying_turbulence_256_strong_surface_waves_statistics.jld2"
filename_weak_waves   = "decaying_turbulence_256_weak_surface_waves_statistics.jld2"
filename_isotropic    = "decaying_turbulence_256_isotropic_statistics.jld2"
filename_rotating     = "decaying_turbulence_256_rotating_statistics.jld2"

emt = FieldTimeSeries(filename_waves,        "e")
est = FieldTimeSeries(filename_strong_waves, "e")
ewt = FieldTimeSeries(filename_weak_waves,   "e")
eit = FieldTimeSeries(filename_isotropic,    "e")
ert = FieldTimeSeries(filename_rotating,     "e")

t = ewt.times
Nt = length(t)

ei = eit[:]
er = ert[:]
ew = ewt[:]
es = est[:]
em = emt[:]

Δt = t[2] - t[1]
ϵi = [-(ei[n+1] - ei[n]) / Δt for n = 1:Nt-1]
ϵr = [-(er[n+1] - er[n]) / Δt for n = 1:Nt-1]
ϵw = [-(ew[n+1] - ew[n]) / Δt for n = 1:Nt-1]
ϵs = [-(es[n+1] - es[n]) / Δt for n = 1:Nt-1]
ϵm = [-(em[n+1] - em[n]) / Δt for n = 1:Nt-1]

ei2 = (ei[1:end-1] .+ ei[2:end]) / 2 
er2 = (er[1:end-1] .+ er[2:end]) / 2
ew2 = (ew[1:end-1] .+ ew[2:end]) / 2
es2 = (es[1:end-1] .+ es[2:end]) / 2
em2 = (em[1:end-1] .+ em[2:end]) / 2

function energy_dissipation(du, u, p, t)
    E = u[1]
    ϵ = u[2]
    Ci = p.Ci
    Cs = p.Cs
    S = p.S

    du[1] = - ϵ
    du[2] = - Ci * ϵ^2 / E - Cs * S * ϵ
    
    return nothing
end

#=
# Simulation
function simulation!(Ci, Cs, S; ti=t[2], tf=t[end])
    E0 = ei[2]
    ϵ0 = (ϵi[1] + ϵi[2]) / 2
    problem = ODEProblem(energy_dissipation, [E0, ϵ0], (ti, tf), (; Ci, Cs, S))
    solution = solve(problem, Tsit5())
    return solution.t, solution.u[1], solution.u[2]
end

Ci = 4.0
Cs = 1.0
ti, ei_s, ϵi_s = simulation!(Ci, Cs, 0)
tw, ew_s, ϵw_s = simulation!(Ci, Cs, 0.01)
tm, em_s, ϵm_s = simulation!(Ci, Cs, 0.05)
ts, es_s, ϵs_s = simulation!(Ci, Cs, 0.1)
tr, er_s, ϵr_s = simulation!(Ci, Cs, 0.1)
=#

function analytical_solution(t, t₀, e₀, ϵ₀, Ci, Cs, S)
    n = abs(1 / (1 - Ci))
    return e₀ * (1 + ϵ₀ / (n * e₀ * Cs * S) * (1 - exp(-Cs * S * (t - t₀))))^(-n)
end

Ci = 5.0
Cs = 1e-4
e₀ = ei[2]
ϵ₀ = 1e-8 #(ϵi[1] + ϵi[2]) / 2
ts = t[2:end]
#ei_s = analytical_solution.(ts, e₀, ϵ₀, Ci, Cs, 0)
ew_s = analytical_solution.(ts, ts[1], e₀, ϵ₀, Ci, Cs, 0.01)
em_s = analytical_solution.(ts, ts[1], e₀, ϵ₀, Ci, Cs, 0.05)
es_s = analytical_solution.(ts, ts[1], e₀, ϵ₀, Ci, Cs, 0.1)
er_s = analytical_solution.(ts, ts[1], e₀, ϵ₀, Ci, Cs, 0.1)
                                          
# For plotting
ei = ei[2:end]
er = er[2:end]
ew = ew[2:end]
es = es[2:end]
em = em[2:end]

ωi = ωi[2:end]
ωr = ωr[2:end]
ωw = ωw[2:end]
ωs = ωs[2:end]
ωm = ωm[2:end]

t2 = (t[1:end-1] .+ t[2:end]) / 2
t = t[2:end]

set_theme!(Theme(fontsize=24))
fig = Figure(resolution=(1200, 500))
ax1 = Axis(fig[1, 1], xscale=log10, yscale=log10, xlabel="Time", ylabel="Volume-integrated kinetic energy")

colors = Makie.wong_colors()
lines!(ax1, t, ei, linewidth=3, color=colors[1], label="Isotropic turbulence")
lines!(ax1, t, ew, linewidth=3, color=colors[2], label="Turbulence beneath weak surface waves")
lines!(ax1, t, em, linewidth=3, color=colors[3], label="Turbulence beneath medium surface waves")
lines!(ax1, t, es, linewidth=3, color=colors[4], label="Turbulence beneath strong surface waves")
lines!(ax1, t, er, linewidth=3, color=colors[5], label="Rotating turbulence")

#lines!(ax1, ts, ei_s, linewidth=1, linestyle=:dash, color=colors[1])
lines!(ax1, ts, ew_s, linewidth=1, linestyle=:dash, color=colors[2])
lines!(ax1, ts, em_s, linewidth=1, linestyle=:dash, color=colors[3])
lines!(ax1, ts, es_s, linewidth=1, linestyle=:dash, color=colors[4])
lines!(ax1, ts, er_s, linewidth=1, linestyle=:dash, color=colors[5])

axislegend(ax1, position=:lb)

display(fig)

# save("decaying_kinetic_energy.png", fig)

