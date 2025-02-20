using GLMakie
using Oceananigans

N = 384

cases = [
    "isotropic",
    "rotating",
    #"strong_surface_waves",
    "medium_surface_waves",
    "weak_surface_waves",
    "very_weak_surface_waves",
]

labels = [
    "Isotropic",
    "Rotating with f = 1/4",
    #"strong_surface_waves",
    "Surface-wave-modulated with ∂z uˢ = z / 2",
    "Surface-wave-modulated with ∂z uˢ = z / 4",
    "Surface-wave-modulated with ∂z uˢ = z / 8",
]

background_vorticities = [
    0.0,
    0.25,
    #1.0,
    0.5,
    0.25,
    0.125,
]

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
fig = Figure(size=(900, 450))
ax1 = Axis(fig[1, 1], xscale=log10, yscale=log10, xlabel="Time", ylabel="Normalized kinetic energy, k / k(t=0)")
#ax2 = Axis(fig[2, 1], xscale=log10, yscale=identity, xlabel="Time", ylabel="Data / model ratio")

for (et, label, Ω) in zip(es, labels, background_vorticities)
    e = view(et.data, 1, 1, 1, :)
    t = et.times    
    Nt = length(t)
    t̃ = view(t, 2:Nt)

    @show e₀ = e[1]
    ẽ = view(e, 2:Nt)
    ln = lines!(ax1, t̃, ẽ ./ e₀, alpha=0.6; linewidth=4, label)

    # Estimate ϵ
    #α = 11/6
    α = 1.76 #11/6
    #α = 1.7 #11/6
    ϵ₀, k₀, t₀ = estimate_dissipation(1, 2, e, t, α)

    if label == "Isotropic"
        e★ = e_isotropic.(t̃, k₀, ϵ₀, t₀, α)
    else
        βΩ = estimate_vorticity_parameter(e, ϵ₀, α)
        @show β = βΩ / Ω
        if !occursin("Rotating", label)
            #βΩ = 0.013868141255849906 * Ω
            #βΩ = 0.0139 * Ω
            βΩ = 0.023 * Ω
        end
        e★ = e_vortical.(t̃, k₀, ϵ₀, t₀, α, βΩ)
    end

    lines!(ax1, t̃, e★ ./ e₀; linewidth=2, linestyle=:dash, color=ln.color)
    #lines!(ax2, t̃, ẽ ./ e★, alpha=0.6; label, linewidth=4)
end

axislegend(ax1, position=:lb)

#ax3 = Axis(fig[3, 1], xscale=log10, yscale=identity, xlabel="Time", ylabel="rms vorticity")

xlims!(ax1, 1e-1, 2e4)
#ylims!(ax2, 0.5, 1.2)

hidespines!(ax1, :t, :r)

display(fig)

save("kinetic_energy.png", fig)

