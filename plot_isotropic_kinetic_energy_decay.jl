using GLMakie
using Oceananigans

Ns = [128, 256]

filenames = ["decaying_turbulence_$(N)_9_isotropic_statistics.jld2" for N in Ns]
es = [FieldTimeSeries(name, "e") for name in filenames]
t = first(es).times

t1 = es[1].times
t2 = es[2].times

function e_isotropic(t, k₀, ϵ₀, t₀, α=11/6)
    γ = 1 / (α - 1)
    if t < t₀
        return NaN
    else
        return k₀ / (1 + ϵ₀ / (γ * k₀) * (t - t₀))^γ
    end
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

set_theme!(Theme(fontsize=24))
fig = Figure(resolution=(1400, 1000))
ax1 = Axis(fig[1, 1], xscale=log10, yscale=log10, xlabel="Time", ylabel="Kinetic energy")
ax2 = Axis(fig[2, 1], xscale=log10, yscale=identity, xlabel="Time", ylabel="Kinetic energy")

for (et, N) in zip(es, Ns)
    e = view(et.data, 1, 1, 1, :)
    t = et.times    
    Nt = length(t)
    t̃ = view(t, 2:Nt)

    e₀ = e[1]
    ẽ = view(e, 2:Nt)
    scatter!(ax1, t̃, ẽ ./ e₀, alpha=0.6, label="data, N = $N")

    # Estimate ϵ
    α = 1.83333
    ϵ₀, k₀, t₀ = estimate_dissipation(1, 100, e, t, α)
    e★ = e_isotropic.(t̃, k₀, ϵ₀, t₀, α)

    lines!(ax1, t̃, e★ ./ e₀, label="α = 11/6, N = $N")
    scatter!(ax2, t̃, ẽ ./ e★, alpha=0.6, label="α = 11/6, N = $N")

    #=
    α = 1.75
    ϵ₀, k₀, t₀ = estimate_dissipation(1, 2, e, t, α)
    e★ = e_isotropic.(t̃, k₀, ϵ₀, t₀, α)

      lines!(ax1, t̃, e★ ./ e₀, label="α = 1.75, N = $N")
    scatter!(ax2, t̃, ẽ  ./ e★, label="α = 1.75, N = $N")

    α = 1.7
    ϵ₀, k₀, t₀ = estimate_dissipation(1, 2, e, t, α)
    e★ = e_isotropic.(t̃, k₀, ϵ₀, t₀, α)

      lines!(ax1, t̃, e★ ./ e₀, label="α = 17/10, N = $N")
    scatter!(ax2, t̃, ẽ  ./ e★, label="α = 17/10, N = $N")
    =#
end

axislegend(ax2, position=:lt)
display(fig)

