#####
##### 2D (kx, kz) energy spectrum of xz slices from decaying turbulence runs.
#####
##### The expected feature: a "barbell" anisotropy in (kx, kz) space, accumulating
##### energy near kx = 0 (depth-alternating jets), in analogy with β-plane turbulence
##### where energy accumulates near ky = 0 (zonal jets). See WC25 §2.2.
#####
##### Numerical hazard: z is Bounded (free-slip), so a naive 2D FFT introduces ringing
##### along the kz axis from the wall discontinuity. We mitigate by (a) subtracting the
##### horizontal mean profile ⟨u⟩(z) and (b) applying a Tukey window in z before FFT.
##### The window scales each Fourier-mode amplitude by a known factor, but for relative
##### comparisons across kx,kz at the same time the bias is uniform.
#####
##### Usage:
#####   julia --project=. plot_xz_spectrum.jl <prefix> [time]
##### where <prefix> is e.g. "decaying_turbulence_256_9_medium_surface_waves"
##### and the script reads <prefix>_xz.jld2.
#####

using Oceananigans
using FFTW
using Statistics
using CairoMakie
using Printf
using JLD2

function tukey_window(N; α=0.25)
    w = ones(N)
    edge = round(Int, α * (N-1) / 2)
    for n = 0:edge-1
        w[n+1] = 0.5 * (1 - cos(π * n / edge))
    end
    for n = 0:edge-1
        w[N-n] = 0.5 * (1 - cos(π * n / edge))
    end
    return w
end

"""
    xz_spectrum(u; window=:tukey, α=0.25, detrend_z=true)

Compute the 2D Fourier energy spectrum |û(kx, kz)|² / 2 of an Nx×Nz field u sampled in
the x-z plane. x is periodic; z is bounded — a window in z is required to avoid wall
ringing. Returns (kx_pos, kz_pos, E2D) where E2D[i, j] is energy density at (kx_pos[i],
kz_pos[j]) summed over conjugate pairs (so it integrates to total energy).
"""
function xz_spectrum(u; window=:tukey, α=0.25, detrend_z=false)
    Nx, Nz = size(u)

    # CAVEAT: do NOT subtract horizontal-mean U(z) — that's exactly the depth-alternating
    # jet structure we want to see in the kx=0 column. The kz-axis "ringing" from
    # bounded-z is suppressed via a Tukey window in z instead.
    u_work = detrend_z ? u .- mean(u, dims=1) : copy(u)

    if window == :tukey
        wz = tukey_window(Nz; α)
        u_work = u_work .* reshape(wz, 1, Nz)
    end

    û = fft(u_work) ./ (Nx * Nz)
    e = abs.(û).^2 ./ 2

    # Shift so kx, kz ∈ [-N/2, N/2-1] order
    e_shift = fftshift(e)
    kx = (-Nx÷2 : Nx÷2-1) .* (2π)
    kz = (-Nz÷2 : Nz÷2-1) .* (2π)

    return kx, kz, e_shift
end

function load_xz_slice(prefix, fieldname, time_target)
    filename = prefix * "_xz.jld2"
    fts = FieldTimeSeries(filename, fieldname)
    t   = fts.times
    n   = argmin(abs.(t .- time_target))
    f3  = Array(interior(fts[n]))
    # 3D slice (indices=(:,1,:)) and Flat-y both produce (Nx, 1, Nz) — drop dim 2.
    return dropdims(f3; dims=2), t[n]
end

function plot_one(prefix, time_target; outname=nothing)
    u_xz, tact = load_xz_slice(prefix, "u", time_target)
    v_xz, _    = load_xz_slice(prefix, "v", time_target)
    w_xz, _    = load_xz_slice(prefix, "w", time_target)

    # Trim staggered fields to a common (Nx, Nz) — drop the wall face row of w (free slip → 0 anyway)
    Nx_min = minimum(size(u, 1) for u in (u_xz, v_xz, w_xz))
    Nz_min = minimum(size(u, 2) for u in (u_xz, v_xz, w_xz))
    u_xz = u_xz[1:Nx_min, 1:Nz_min]
    v_xz = v_xz[1:Nx_min, 1:Nz_min]
    w_xz = w_xz[1:Nx_min, 1:Nz_min]

    kx, kz, Eu = xz_spectrum(u_xz)
    _,  _,  Ev = xz_spectrum(v_xz)
    _,  _,  Ew = xz_spectrum(w_xz)
    E = Eu .+ Ev .+ Ew  # full KE spectrum; v carries the cross-wave momentum (Fig 3 of WC25)

    set_theme!(Theme(fontsize=20))
    fig = Figure(size=(900, 500))

    # Real-space u panel
    ax1 = Axis(fig[1, 1], aspect=1, title=@sprintf("u(x, z) at t = %.0f", tact),
               xlabel="x", ylabel="z")
    ulim = maximum(abs, u_xz)
    Nx, Nz = size(u_xz)
    x = range(0, 1, length=Nx); z = range(0, 1, length=Nz)
    heatmap!(ax1, x, z, u_xz; colormap=:balance, colorrange=(-ulim, ulim))

    # Spectrum: log-log, restrict to positive kx and full kz, show "barbell" anisotropy.
    Eplot = max.(E, 1e-30)
    logE = log10.(Eplot)
    vmin = maximum(logE) - 5
    vmax = maximum(logE)

    ax2 = Axis(fig[1, 2], aspect=1, title="log₁₀ E(kx, kz)",
               xlabel="kx / 2π", ylabel="kz / 2π")
    hm  = heatmap!(ax2, kx ./ (2π), kz ./ (2π), logE;
                   colormap=:viridis, colorrange=(vmin, vmax))
    Colorbar(fig[1, 3], hm)

    Klim = 32   # show out to the IC peak wavenumber k₀/2π = 32
    xlims!(ax2, -Klim, Klim); ylims!(ax2, -Klim, Klim)

    # Overlay the Rhines-like wave-arrest scale |kz| ~ √(∂²u^S / U) using
    # ∂²u^S = 0.5 (medium waves) and a representative U from the field.
    if occursin("medium", prefix) || occursin("xz25", prefix)
        ∂²uˢ = 0.5
        U    = sqrt(2 * mean(u_xz.^2 .+ w_xz.^2))  # rms speed proxy
        if U > 0
            kR = sqrt(∂²uˢ / U) / (2π)
            lines!(ax2, [-Klim, Klim], [ kR,  kR]; color=:red, linestyle=:dash)
            lines!(ax2, [-Klim, Klim], [-kR, -kR]; color=:red, linestyle=:dash)
        end
    end

    name = something(outname, prefix * @sprintf("_xz_spectrum_t%05d.png", round(Int, tact)))
    save(name, fig)
    @info "wrote $name"
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    prefix = length(ARGS) >= 1 ? ARGS[1] : "decaying_turbulence_256_9_medium_surface_waves"
    tval   = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1e4
    plot_one(prefix, tval)
end
