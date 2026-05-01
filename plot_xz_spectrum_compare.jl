#####
##### Side-by-side (kx, kz) energy spectra for the medium-waves and isotropic cases,
##### with 1D marginal spectra to make the barbell anisotropy visible.
#####
##### Usage:
#####   julia --project=. plot_xz_spectrum_compare.jl [time]
#####
##### Sources data from xz-slice files. If a *_fields.jld2 file is present at the
##### requested time, it is preferred (full 3D — we y-average the 2D FFT for a much
##### smoother spectrum).
#####

include("plot_xz_spectrum.jl")

using CairoMakie

"""
    spectrum_from_3d_field(prefix, time_target)

Compute (kx, kz) spectrum from a full 3D field by 2D-FFT-ing each y-slice and averaging.
Falls back to the single-y xz slice if the fields file isn't available at this time.
"""
function spectrum_from_3d_field(prefix, time_target; α=0.25)
    fields_file = prefix * "_fields.jld2"
    if isfile(fields_file)
        try
            ut = FieldTimeSeries(fields_file, "u")
            vt = FieldTimeSeries(fields_file, "v")
            wt = FieldTimeSeries(fields_file, "w")
            t  = ut.times
            n  = argmin(abs.(t .- time_target))
            if abs(t[n] - time_target) < 1.0
                u = Array(interior(ut[n]))
                v = Array(interior(vt[n]))
                w = Array(interior(wt[n]))
                Nx_min = minimum(size.((u, v, w), 1))
                Ny_min = minimum(size.((u, v, w), 2))
                Nz_min = minimum(size.((u, v, w), 3))
                u = u[1:Nx_min, 1:Ny_min, 1:Nz_min]
                v = v[1:Nx_min, 1:Ny_min, 1:Nz_min]
                w = w[1:Nx_min, 1:Ny_min, 1:Nz_min]

                E_avg = nothing
                kx = kz = nothing
                for j = 1:Ny_min
                    kx, kz, Eu = xz_spectrum(@view u[:, j, :]; α)
                    _,  _,  Ev = xz_spectrum(@view v[:, j, :]; α)
                    _,  _,  Ew = xz_spectrum(@view w[:, j, :]; α)
                    E = Eu .+ Ev .+ Ew
                    E_avg = isnothing(E_avg) ? E : E_avg .+ E
                end
                return kx, kz, E_avg ./ Ny_min, t[n], "3D (Ny=$Ny_min y-slices avg)"
            end
        catch e
            @warn "Couldn't load 3D fields, falling back to xz slice: $e"
        end
    end
    # Fallback to single y=1 xz slice
    u_xz, tact = load_xz_slice(prefix, "u", time_target)
    v_xz, _    = load_xz_slice(prefix, "v", time_target)
    w_xz, _    = load_xz_slice(prefix, "w", time_target)
    Nx_min = minimum(size(u, 1) for u in (u_xz, v_xz, w_xz))
    Nz_min = minimum(size(u, 2) for u in (u_xz, v_xz, w_xz))
    u_xz = u_xz[1:Nx_min, 1:Nz_min]
    v_xz = v_xz[1:Nx_min, 1:Nz_min]
    w_xz = w_xz[1:Nx_min, 1:Nz_min]
    kx, kz, Eu = xz_spectrum(u_xz; α)
    _,  _,  Ev = xz_spectrum(v_xz; α)
    _,  _,  Ew = xz_spectrum(w_xz; α)
    return kx, kz, Eu .+ Ev .+ Ew, tact, "xz slice (y=1)"
end

function plot_compare(prefix_iso, prefix_med, time_target; outname="xz_spectrum_compare.png", Klim=24)
    kx_i, kz_i, E_i, t_i, source_i = spectrum_from_3d_field(prefix_iso, time_target)
    kx_m, kz_m, E_m, t_m, source_m = spectrum_from_3d_field(prefix_med, time_target)

    @info "isotropic source: $source_i, t=$t_i"
    @info "medium    source: $source_m, t=$t_m"

    # 1D marginal spectra (sum over the OTHER dimension)
    # E(kx) = Σ_kz E(kx, kz);  E(kz) = Σ_kx E(kx, kz)
    Ei_kx = vec(sum(E_i, dims=2))
    Ei_kz = vec(sum(E_i, dims=1))
    Em_kx = vec(sum(E_m, dims=2))
    Em_kz = vec(sum(E_m, dims=1))

    # Take only positive wavenumbers (fold conjugate pair)
    Nx = length(kx_i); Nz = length(kz_i)
    pos_kx_idx = (Nx ÷ 2 + 1):Nx          # kx=0 then positive
    pos_kz_idx = (Nz ÷ 2 + 1):Nz
    kx_pos = kx_i[pos_kx_idx]
    kz_pos = kz_i[pos_kz_idx]
    Ei_kx_pos = Ei_kx[pos_kx_idx]; Ei_kz_pos = Ei_kz[pos_kz_idx]
    Em_kx_pos = Em_kx[pos_kx_idx]; Em_kz_pos = Em_kz[pos_kz_idx]

    set_theme!(Theme(fontsize=22))
    fig = Figure(size=(1400, 1500))

    # Common colorbar range from both
    Eall = vcat(vec(E_i), vec(E_m))
    vmax = log10(maximum(Eall))
    vmin = vmax - 5

    ax_iso = Axis(fig[1, 1], aspect=1, title=@sprintf("isotropic   t=%.0f", t_i),
                  xlabel="kx / 2π", ylabel="kz / 2π")
    hm_iso = heatmap!(ax_iso, kx_i ./ (2π), kz_i ./ (2π), log10.(max.(E_i, 1e-30));
                      colormap=:viridis, colorrange=(vmin, vmax))
    xlims!(ax_iso, -Klim, Klim); ylims!(ax_iso, -Klim, Klim)

    ax_med = Axis(fig[1, 2], aspect=1, title=@sprintf("medium waves   t=%.0f", t_m),
                  xlabel="kx / 2π", ylabel="kz / 2π")
    heatmap!(ax_med, kx_m ./ (2π), kz_m ./ (2π), log10.(max.(E_m, 1e-30));
             colormap=:viridis, colorrange=(vmin, vmax))
    xlims!(ax_med, -Klim, Klim); ylims!(ax_med, -Klim, Klim)

    Colorbar(fig[1, 3], hm_iso; label="log₁₀ E(kx, kz)")

    # Rhines-like overlay on medium panel
    ∂²uˢ = 0.5
    U    = sqrt(2 * sum(E_m))     # rms speed proxy
    if U > 0
        kR = sqrt(∂²uˢ / U) / (2π)
        lines!(ax_med, [-Klim, Klim], [ kR,  kR]; color=:red, linestyle=:dash, linewidth=2)
        lines!(ax_med, [-Klim, Klim], [-kR, -kR]; color=:red, linestyle=:dash, linewidth=2)
        text!(ax_med, -Klim + 1, kR + 0.7;
              text=@sprintf("kᴿ ≈ %.1f", kR),
              color=:red, fontsize=18)
    end

    # 1D marginals: kx (left) and kz (right) cuts
    ax_kx = Axis(fig[2, 1:2], xscale=log10, yscale=log10,
                 title="E(kx) = Σ_kz E(kx, kz)",
                 xlabel="kx / 2π", ylabel="E(kx)")
    lines!(ax_kx, kx_pos[2:end] ./ (2π), Ei_kx_pos[2:end]; label="isotropic", color=:gray)
    lines!(ax_kx, kx_pos[2:end] ./ (2π), Em_kx_pos[2:end]; label="medium waves", color=:steelblue, linewidth=2)
    axislegend(ax_kx; position=:lb)

    ax_kz = Axis(fig[3, 1:2], xscale=log10, yscale=log10,
                 title="E(kz) = Σ_kx E(kx, kz)   ←  energy in (depth-alternating jets) for medium waves",
                 xlabel="kz / 2π", ylabel="E(kz)")
    lines!(ax_kz, kz_pos[2:end] ./ (2π), Ei_kz_pos[2:end]; label="isotropic", color=:gray)
    lines!(ax_kz, kz_pos[2:end] ./ (2π), Em_kz_pos[2:end]; label="medium waves", color=:steelblue, linewidth=2)
    axislegend(ax_kz; position=:lb)

    rowsize!(fig.layout, 1, Relative(0.5))
    save(outname, fig, px_per_unit=2)
    @info "wrote $outname"
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    tval   = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1e4
    pi_iso = "decaying_turbulence_256_9_isotropic"
    pi_med = "decaying_turbulence_256_9_medium_surface_waves"
    plot_compare(pi_iso, pi_med, tval; outname=@sprintf("xz_spectrum_compare_t%05d.png", round(Int, tval)))
end
