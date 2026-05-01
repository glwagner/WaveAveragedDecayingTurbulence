#####
##### Sanity-check the (kx, kz) spectrum recipe on a synthetic field that
##### intentionally has a "barbell" anisotropy and a wall-bounded mean profile.
##### If the plotter is right, the synthetic barbell should show up cleanly along
##### the kz axis, with no kz-axis ringing.
#####

include("../plot_xz_spectrum.jl")

using FFTW, Statistics, CairoMakie, Printf

Nx = Nz = 256
x = range(0, 1, length=Nx)
z = range(0, 1, length=Nz)

# Build a field with energy concentrated at low kx (jets along z, like depth-alternating jets).
# u(x, z) = sum over a few jet modes that are k_x≈0, k_z = m·2π for m = ±1..±6
import Random; Random.seed!(7)
u = zeros(Nx, Nz)
for m = 1:6
    amp = 1.0 / m
    u .+= amp .* cos.(m .* 2π .* z)' .* (1 .+ 0.05 .* cos.(2π .* x))  # tiny x dependence
end

# Sprinkle a small isotropic background and a wall-anchored mean profile (mimics imperfect
# subtraction of a boundary layer).
for k_iso = 1:8, l_iso = 1:8
    u .+= 0.02 .* sin.(2π * k_iso .* x .+ rand()) .* sin.(2π * l_iso .* z .+ rand())'
end
mean_profile = 0.1 .* (z .- 0.5)  # linear wall-to-wall trend, will dominate kx=0 column if not removed
u .+= mean_profile'

# 1. With detrending + Tukey window (the recipe).
kx, kz, E_clean = xz_spectrum(u; window=:tukey, α=0.25, detrend_z=true)

# 2. Without detrending or windowing (the bad recipe — should show kz-axis artifact).
kx, kz, E_bad   = xz_spectrum(u; window=:none,  α=0.0,  detrend_z=false)

set_theme!(Theme(fontsize=18))
fig = Figure(size=(1200, 500))

ax0 = Axis(fig[1, 1], aspect=1, title="synthetic u(x, z)", xlabel="x", ylabel="z")
heatmap!(ax0, x, z, u; colormap=:balance, colorrange=(-2, 2))

vmax = log10(max(maximum(E_clean), maximum(E_bad)))
vmin = vmax - 6

ax1 = Axis(fig[1, 2], aspect=1, title="bad recipe — naive 2D FFT",
           xlabel="kx / 2π", ylabel="kz / 2π")
hm1 = heatmap!(ax1, kx ./ (2π), kz ./ (2π), log10.(max.(E_bad, 1e-30));
               colormap=:viridis, colorrange=(vmin, vmax))
xlims!(ax1, -16, 16); ylims!(ax1, -16, 16)

ax2 = Axis(fig[1, 3], aspect=1, title="recipe (detrend + Tukey)",
           xlabel="kx / 2π", ylabel="kz / 2π")
hm2 = heatmap!(ax2, kx ./ (2π), kz ./ (2π), log10.(max.(E_clean, 1e-30));
               colormap=:viridis, colorrange=(vmin, vmax))
xlims!(ax2, -16, 16); ylims!(ax2, -16, 16)
Colorbar(fig[1, 4], hm2; label="log₁₀ E")

save("synthetic_spectrum_check.png", fig)
@info "wrote synthetic_spectrum_check.png"
