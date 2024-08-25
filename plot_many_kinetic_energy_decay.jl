using GLMakie
using Oceananigans

N = 256

cases = ["isotropic",
         "rotating",
         #"strong_surface_waves",
         "medium_surface_waves",
         "surface_waves"]

labels = ["isotropic",
          "rotating",
         # "strong_surface_waves",
          "medium_surface_waves",
          "weak_surface_waves"]

filenames = ["decaying_turbulence_$(N)_9_$(case)_statistics.jld2" for case in cases]
es = [FieldTimeSeries(name, "e") for name in filenames]

function e_isotropic(t, k₀, ϵ₀, t₀, α)
    γ = 1 / (α - 1)
    if t < t₀
        return NaN
    else
        return k₀ / (1 + ϵ₀ / (γ * k₀) * (t - t₀))^γ
    end
end

function e_vortical(t, k₀, ϵ₀, t₀, α, βΩ)
    γ = 1 / (α - 1)
    return k₀ / (1 + ϵ₀ / (γ * βΩ * k₀) * (1 - exp(-βΩ * (t - t₀))))^γ
end

function estimate_dissipation(n₀, n₊, e, t, α)
    γ = 1 / (α - 1)

    k₀ = e[n₀]
    t₀ = t[n₀]
    k₊ = e[n₊]
    t₊ = t[n₊]
    Δt = (t₊ - t₀)

    r = k₀ / k₊
    ϵ₀ = γ * k₀ / Δt * (r^(1/γ) - 1)

    return ϵ₀, k₀, t₀
end

function estimate_vorticity_parameter(e, ϵ₀, α)
    γ = 1 / (α - 1)

    k₀ = e[1]
    k∞ = e[end]

    σ = k₀ / k∞ 

    βΩ = ϵ₀ / (γ * k₀ * (σ^(1/γ) - 1))

    @show βΩ

    return βΩ
end

set_theme!(Theme(fontsize=24))
fig = Figure(size=(1400, 1000))
ax1 = Axis(fig[1, 1], xscale=log10, yscale=log10, xlabel="Time", ylabel="Kinetic energy")
ax2 = Axis(fig[2, 1], xscale=log10, yscale=identity, xlabel="Time", ylabel="Data / model ratio")

for (et, label) in zip(es, labels)
    e = view(et.data, 1, 1, 1, :)
    t = et.times    
    Nt = length(t)
    t̃ = view(t, 2:Nt)

    e₀ = e[1]
    ẽ = view(e, 2:Nt)
    ln = lines!(ax1, t̃, ẽ ./ e₀, alpha=0.6; linewidth=6, label)

    # Estimate ϵ
    α = 11/6
    ϵ₀, k₀, t₀ = estimate_dissipation(1, 100, e, t, α)

    if label == "isotropic"
        e★ = e_isotropic.(t̃, k₀, ϵ₀, t₀, α)
    else
        βΩ = estimate_vorticity_parameter(e, ϵ₀, α)
        e★ = e_vortical.(t̃, k₀, ϵ₀, t₀, α, βΩ)
    end

    lines!(ax1, t̃, e★ ./ e₀; label, linewidth=2, linestyle=:dash, color=ln.color)
    lines!(ax2, t̃, ẽ ./ e★, alpha=0.6; label, linewidth=4)
end

axislegend(ax2, position=:lt)

ax3 = Axis(fig[3, 1], xscale=log10, yscale=identity, xlabel="Time", ylabel="rms vorticity")

N = 256
case = "rotating"
filename = "decaying_turbulence_$(N)_9_$(case)_statistics.jld2"
ω²t = FieldTimeSeries(filename, "ω²")
t = first(es).times
Ω = sqrt.(ω²t[:])

lines!(ax3, t[2:end], Ω[2:end], linewidth=4)

t0 = 1
t1 = 1e2

#=
xlims!(ax1, t0, t1)
xlims!(ax2, t0, t1)
xlims!(ax3, t0, t1)
=#

ylims!(ax2, 0.9, 1.2)

display(fig)

