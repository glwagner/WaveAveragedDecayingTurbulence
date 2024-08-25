using GLMakie
using Oceananigans

N = 512
filename_isotropic       = "decaying_turbulence_$(N)_isotropic_statistics.jld2"





#filename_waves           = "decaying_turbulence_$(N)_surface_waves_statistics.jld2"
#filename_rotating        = "decaying_turbulence_$(N)_rotating_statistics.jld2"
#filename_strong_waves    = "decaying_turbulence_$(N)_strong_surface_waves_statistics.jld2"
#filename_weak_waves      = "decaying_turbulence_$(N)_weak_surface_waves_statistics.jld2"
#filename_very_weak_waves = "decaying_turbulence_$(N)_very_weak_surface_waves_statistics.jld2"


eit = FieldTimeSeries(filename_isotropic, "e")
ert = FieldTimeSeries(filename_rotating, "e")
est = FieldTimeSeries(filename_strong_waves, "e")
emt = FieldTimeSeries(filename_waves, "e")
#ewt = FieldTimeSeries(filename_weak_waves, "e")
#evt = FieldTimeSeries(filename_very_weak_waves, "e")

n₀ = 1

ei = eit[:]
er = ert[:]
em = emt[:]
es = est[:]
#ew = ewt[:]
#ev = evt[:]

t = eit.times
Nt = length(t)

ei = ei[n₀:end]
er = er[n₀:end]
em = em[n₀:end]
es = es[n₀:end]
  t = t[n₀:end]

#ew = ew[n₀:end]
#ev = ev[n₀:end]

set_theme!(Theme(fontsize=24))
fig = Figure(resolution=(1400, 1000))
ax1 = Axis(fig[1, 1], xscale=log10, yscale=log10, xlabel="Time", ylabel="Kinetic energy")
#ax1 = Axis(fig[1, 1], yscale=log10, xlabel="Time", ylabel="Kinetic energy")

colors = Makie.wong_colors()
scatter!(ax1, t[2:end], ei[2:end] ./ ei[1], linewidth=3, color=colors[1], label="Isotropic turbulence")
scatter!(ax1, t[2:end], em[2:end] ./ em[1], linewidth=3, color=colors[3], label="Turbulence beneath medium surface waves")
scatter!(ax1, t[2:end], er[2:end] ./ er[1], linewidth=3, color=colors[3], label="Rotating turbulence")
scatter!(ax1, t[2:end], es[2:end] ./ es[1], linewidth=3, color=colors[3], label="Turbulence beneath strong surface waves")

#lines!(ax1, t, es, linewidth=3, color=colors[2], label="Turbulence beneath strong surface waves")
#lines!(ax1, t, ew, linewidth=3, color=colors[4], label="Turbulence beneath weak surface waves")
#lines!(ax1, t, ev, linewidth=3, color=colors[5], label="Turbulence beneath very weak surface waves")
#lines!(ax1, t, er, linewidth=3, color=colors[6], label="Rotating turbulence")

ϵ(n) = - (em[n] - em[n-2]) / (t[n] - t[n-2])

α = 1.76 # 11/6

# estimate epsilon
n₀ = 1
k₀ = ei[n₀]
t₀ = t[n₀]

n₁ = 2
k₁ = ei[n₁]
t₁ = t[n₁]
ϵ₀ = k₀ / (t₁ * (1 - α)) * (1 - (k₁ / k₀)^(1 - α))
@show ϵ₀ ϵ(3)

γ = 1 - α
μ = 1 / (1 - α)
em1 = @. k₀ * (α - 1)^μ * (1 + ϵ₀ / k₀ * (t - t₀))^μ
lines!(ax1, t[2:end], em1[2:end] ./ ei[1], linestyle=:dot, linewidth=1, color=(:black, 0.8))

N = 1001
k∞ = er[N]
@show βS = (α - 1) * (ϵ₀ / k₀) / ((k∞/k₀)^γ - 1)
em2 = @. k₀ * (1 - (α - 1) / βS * ϵ₀ / k₀ * (exp(-βS * (t - t₀)) - 1))^μ
lines!(ax1, t[2:end], em2[2:end] ./ ei[1], linestyle=:solid, linewidth=1, color=(:black, 0.8))

N = 150
k∞ = es[N]
@show βS = (α - 1) * (ϵ₀ / k₀) / ((k∞/k₀)^γ - 1)
em2 = @. k₀ * (1 - (α - 1) / βS * ϵ₀ / k₀ * (exp(-βS * (t - t₀)) - 1))^μ
lines!(ax1, t[2:end], em2[2:end] ./ ei[1], linestyle=:solid, linewidth=1, color=(:black, 0.8))

N = 500
k∞ = em[N]
@show βS = (α - 1) * (ϵ₀ / k₀) / ((k∞/k₀)^γ - 1)
em2 = @. k₀ * (1 - (α - 1) / βS * ϵ₀ / k₀ * (exp(-βS * (t - t₀)) - 1))^μ
lines!(ax1, t[2:end], em2[2:end] ./ ei[1], linestyle=:solid, linewidth=1, color=(:black, 0.8))

#=
ei★ = ei[100]
t★ = t[100]
ei★ = ei[1]
t★ = t[1]
e⁺ = @. ei★ * (t / t★)^(-6/5)
e⁻ = @. ei★ * (t / t★)^(-10/7)

# E ~ t⁻ⁿ
# ϵ ~ t⁻ⁿ⁻¹
# ω ~ ϵ / E ~ t⁻¹
ω₀ = ωi[100] 
ω★ = @. ω₀ * (t / t★)^(-4/5)
lines!(ax2, t, ω★, linewidth=6, color=(:black, 0.2))
band!(ax1, t, e⁺, e⁻, color=(:black, 0.2))
lines!(ax1, t, e⁺, linestyle=:dash, linewidth=1, color=(:black, 0.8))
lines!(ax1, t, e⁻, linestyle=:solid, linewidth=1, color=(:purple, 0.8))
=#

axislegend(ax1, position=:lb)
display(fig)

# save("decaying_kinetic_energy_$(cnt).png", fig)

#=
filename_waves        = "decaying_turbulence_256_surface_waves_statistics.jld2"
filename_strong_waves = "decaying_turbulence_256_strong_surface_waves_statistics.jld2"
filename_weak_waves   = "decaying_turbulence_256_weak_surface_waves_statistics.jld2"
filename_rotating     = "decaying_turbulence_256_rotating_statistics.jld2"

emt = FieldTimeSeries(filename_waves,        "e")
est = FieldTimeSeries(filename_strong_waves, "e")
ewt = FieldTimeSeries(filename_weak_waves,   "e")
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
=#
                                          
# For plotting
#er = er[2:end]
#ew = ew[2:end]
#es = es[2:end]
#em = em[2:end]

#t2 = (t[1:end-1] .+ t[2:end]) / 2

#lines!(ax1, t, ew, linewidth=3, color=colors[2], label="Turbulence beneath weak surface waves")
#lines!(ax1, t, em, linewidth=3, color=colors[3], label="Turbulence beneath medium surface waves")
#lines!(ax1, t, es, linewidth=3, color=colors[4], label="Turbulence beneath strong surface waves")
#lines!(ax1, t, er, linewidth=3, color=colors[5], label="Rotating turbulence")

#=
lines!(ax2, t2, max.(1e-16, ϵi), linewidth=3, color=colors[1], label="Isotropic turbulence")
#lines!(ax2, t2, max.(1e-16, ϵw), linewidth=3, color=colors[2], label="Turbulence beneath weak surface waves")
#lines!(ax2, t2, max.(1e-16, ϵm), linewidth=3, color=colors[3], label="Turbulence beneath medium surface waves")
#lines!(ax2, t2, max.(1e-16, ϵs), linewidth=3, color=colors[4], label="Turbulence beneath strong surface waves")
lines!(ax2, t2, max.(1e-16, ϵr), linewidth=3, color=colors[5], label="Rotating turbulence")

#lines!(ax1, ts, ei_s, linewidth=1, linestyle=:dash, color=colors[1])
lines!(ax1, ts, ew_s, linewidth=1, linestyle=:dash, color=colors[2])
lines!(ax1, ts, em_s, linewidth=1, linestyle=:dash, color=colors[3])
lines!(ax1, ts, es_s, linewidth=1, linestyle=:dash, color=colors[4])
lines!(ax1, ts, er_s, linewidth=1, linestyle=:dash, color=colors[5])
=#

# save("decaying_kinetic_energy.png", fig)

