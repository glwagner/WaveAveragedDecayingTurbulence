using FourierFlows
using FFTW
using Oceananigans
using Statistics
using GLMakie

filename = "decaying_turbulence_512_isotropic_fields.jld2"

# times = [0, 100, 1000, 10000]
# ut = FieldTimeSeries(filename, "u", times=[0, 100, 1000])

spectral_grid = ThreeDGrid(nx=512, Lx=1)
k = spectral_grid.kr

En = []

for n = [1, 2, 3]
    u = interior(ut[n])
    û = rfft(u, 1)

    e = @. abs(û)^2 / 2
    E = mean(e, dims=(2, 3))

    push!(En, E)
end

# Intended initial condition
k₀ = 16 * 2π
k′ = @. k / k₀
uᵢ = @. k′ * sqrt(exp(2 - (k′ + 1)^2 / 2))
Eᵢ = @. uᵢ^2 / 2
Eᵢ .*= E[2, 1, 1] / Eᵢ[2, 1, 1]

fig = Figure()
ax = Axis(fig[1, 1], xscale=log10, yscale=log10)
ylims!(ax, 1e-30, 1e5)

lines!(ax, k[2:end, 1, 1], En[1][2:end, 1, 1])
lines!(ax, k[2:end, 1, 1], En[2][2:end, 1, 1])
lines!(ax, k[2:end, 1, 1], En[3][2:end, 1, 1])

#lines!(ax, k[2:end, 1, 1], Eᵢ[2:end, 1, 1])

display(fig)

# ⟨ u(t) u(t + τ) ⟩
