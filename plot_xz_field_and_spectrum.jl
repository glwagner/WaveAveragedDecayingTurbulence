#####
##### Real-space xz vorticity (η) and (kx, kz) energy spectrum side-by-side for
##### the isotropic and medium-waves cases, plus a wave/isotropic anisotropy panel
##### that exposes the wave-arrest "barbell" shape directly.
#####
##### Layout (3 rows × 2 cols + colorbar col):
#####   row 1: η(x, z) — isotropic | medium waves      (real-space coherent structures)
#####   row 2: log₁₀ E(kx, kz) — isotropic | medium waves   (raw spectra)
#####   row 3: log₁₀ [E_medium / E_isotropic]
#####          (this is where the barbell is unmistakable: a peanut-shaped
#####           wave-induced enhancement aligned with the kz axis)
#####
##### Usage:
#####   julia --project=. plot_xz_field_and_spectrum.jl [time]
#####

include("plot_xz_spectrum_compare.jl")

using CairoMakie

"""
    eta_xz_at(prefix, t_target)

Mid-y xz slice of η = ∂ₓw − ∂_z u from the 3D fields snapshot at t_target.
"""
function eta_xz_at(prefix, t_target)
    fields_file = prefix * "_fields.jld2"
    if isfile(fields_file)
        ut = FieldTimeSeries(fields_file, "u")
        wt = FieldTimeSeries(fields_file, "w")
        n = argmin(abs.(ut.times .- t_target))
        if abs(ut.times[n] - t_target) < 1.0
            u3 = Array(interior(ut[n])); w3 = Array(interior(wt[n]))
            Ny = size(u3, 2); j = max(1, Ny ÷ 2)
            Nz = min(size(u3, 3), size(w3, 3))
            Nx = min(size(u3, 1), size(w3, 1))
            u_xz = u3[1:Nx, j, 1:Nz]; w_xz = w3[1:Nx, j, 1:Nz]
            Δx = 1 / Nx; Δz = 1 / Nz
            ∂x_w = (circshift(w_xz, (-1, 0)) .- circshift(w_xz, (1, 0))) ./ (2Δx)
            ∂z_u = similar(u_xz)
            ∂z_u[:, 2:end-1] .= (u_xz[:, 3:end] .- u_xz[:, 1:end-2]) ./ (2Δz)
            ∂z_u[:, 1] .= (u_xz[:, 2] .- u_xz[:, 1]) ./ Δz
            ∂z_u[:, end] .= (u_xz[:, end] .- u_xz[:, end-1]) ./ Δz
            return ∂x_w .- ∂z_u, ut.times[n]
        end
    end
    return load_xz_slice(prefix, "η", t_target)
end

function plot_field_and_spectrum(prefix_iso, prefix_med, time_target;
                                 outname=@sprintf("xz_field_and_spectrum_t%05d.png", round(Int, time_target)),
                                 Klim=24, Klim_aniso=8, log_decades=1.5)

    η_i, t_i = eta_xz_at(prefix_iso, time_target)
    η_m, t_m = eta_xz_at(prefix_med, time_target)
    Nx = min(size(η_i,1), size(η_m,1)); Nz = min(size(η_i,2), size(η_m,2))
    η_i = η_i[1:Nx, 1:Nz]; η_m = η_m[1:Nx, 1:Nz]
    x = range(0, 1, length=Nx); z = range(0, 1, length=Nz)

    # Percentile-based color range — common across both panels
    p99 = quantile(vcat(vec(abs.(η_i)), vec(abs.(η_m))), 0.99)
    ηlim = p99
    @info "η colorrange = ±$(round(ηlim, sigdigits=3)) (99th percentile)"

    # Spectra (y-averaged from full 3D fields)
    kx, kz, E_i, _, source_i = spectrum_from_3d_field(prefix_iso, time_target)
    _,  _,  E_m, _, source_m = spectrum_from_3d_field(prefix_med, time_target)
    @info "isotropic spectrum source: $source_i"
    @info "medium    spectrum source: $source_m"

    # Tight log range — the kᵦ peak is only ~1 decade above the kz=0 hollow on
    # the dumbbell, so a wide log range buries it. log_decades=1.5 gives clear
    # contrast across the dumbbell without saturating everything.
    vmax = log10(maximum(vcat(vec(E_i), vec(E_m))))
    vmin = vmax - log_decades
    @info "spectrum log range: [$(round(vmin, digits=2)), $(round(vmax, digits=2))]"

    # Anisotropy: log ratio E_medium / E_isotropic. Floor very-low-energy regions of
    # both spectra (near grid scale) so we don't divide noise by noise.
    floor_E = 1e-6 * max(maximum(E_i), maximum(E_m))
    R = log10.(max.(E_m, floor_E) ./ max.(E_i, floor_E))
    rlim = 1.5

    set_theme!(Theme(fontsize=22))
    fig = Figure(size=(1500, 1700))

    # --- Row 1: real-space η(x, z) ---
    ax_iso_r = Axis(fig[1, 1]; aspect=1, title="isotropic",
                    ylabel="z", xticks=[0, 0.5, 1], yticks=[0, 0.5, 1])
    hm_r = heatmap!(ax_iso_r, x, z, η_i; colormap=:balance, colorrange=(-ηlim, ηlim))
    hidexdecorations!(ax_iso_r, ticks=false)

    ax_med_r = Axis(fig[1, 2]; aspect=1, title="medium surface waves",
                    xticks=[0, 0.5, 1], yticks=[0, 0.5, 1])
    heatmap!(ax_med_r, x, z, η_m; colormap=:balance, colorrange=(-ηlim, ηlim))
    hidexdecorations!(ax_med_r, ticks=false); hideydecorations!(ax_med_r, ticks=false)

    Colorbar(fig[1, 3], hm_r; label="η = ∂ₓw − ∂_z u")

    # --- Row 2: log E(kx, kz), tight range, zoomed to barbell scale ---
    ax_iso_s = Axis(fig[2, 1]; aspect=1, ylabel="kz / 2π")
    hm_s = heatmap!(ax_iso_s, kx ./ (2π), kz ./ (2π), log10.(max.(E_i, 1e-30));
                    colormap=:viridis, colorrange=(vmin, vmax))
    xlims!(ax_iso_s, -Klim_aniso, Klim_aniso); ylims!(ax_iso_s, -Klim_aniso, Klim_aniso)
    hidexdecorations!(ax_iso_s, ticks=false)

    ax_med_s = Axis(fig[2, 2]; aspect=1)
    heatmap!(ax_med_s, kx ./ (2π), kz ./ (2π), log10.(max.(E_m, 1e-30));
             colormap=:viridis, colorrange=(vmin, vmax))
    xlims!(ax_med_s, -Klim_aniso, Klim_aniso); ylims!(ax_med_s, -Klim_aniso, Klim_aniso)
    hidexdecorations!(ax_med_s, ticks=false); hideydecorations!(ax_med_s, ticks=false)

    # Overlay the kᵦ peanut on the medium-waves raw spectrum so the user can
    # see the wave-arrest envelope directly on the energy distribution.
    ∂²uˢ_overlay = 0.5
    U_overlay = sqrt(2 * sum(E_m))
    if U_overlay > 0
        kR_overlay = sqrt(∂²uˢ_overlay / U_overlay) / (2π)
        θ_o = range(0, 2π, length=400)
        kβ_o = @. kR_overlay * sqrt(abs(cos(θ_o)))
        lines!(ax_med_s, kβ_o .* cos.(θ_o), kβ_o .* sin.(θ_o);
               color=:red, linestyle=:dash, linewidth=2)
    end

    Colorbar(fig[2, 3], hm_s; label="log₁₀ E(kx, kz)")

    # --- Row 3: anisotropy log10(E_medium / E_isotropic) — the barbell ---
    ax_aniso = Axis(fig[3, 1:2]; aspect=DataAspect(),
                    xlabel="kx / 2π", ylabel="kz / 2π",
                    title="wave-induced anisotropy:  log₁₀ ( E_medium / E_isotropic )")
    hm_aniso = heatmap!(ax_aniso, kx ./ (2π), kz ./ (2π), R;
                        colormap=:balance, colorrange=(-rlim, rlim))
    xlims!(ax_aniso, -Klim_aniso, Klim_aniso); ylims!(ax_aniso, -Klim_aniso, Klim_aniso)

    # Overlay the wave-arrest scale kᴿ as a peanut/barbell-shaped boundary
    ∂²uˢ = 0.5
    U = sqrt(2 * sum(E_m))
    if U > 0
        kR = sqrt(∂²uˢ / U) / (2π)
        # Rhines / Vallis-Maltrud anisotropy boundary in (kx, kz):
        #     k_β(θ) = kᴿ √|cos θ|,    θ measured from kx axis.
        # → peanut-shaped, lobes along kx, pinched at kz axis (θ=±π/2).
        # Inside this curve, "Rossby"-like waves dominate over turbulence and
        # energy is suppressed; it accumulates outside, especially near kx=0
        # (the depth-alternating-jet axis).
        θ = range(0, 2π, length=400)
        kβ = @. kR * sqrt(abs(cos(θ)))
        kx_b = @. kβ * cos(θ) * (2π)
        kz_b = @. kβ * sin(θ) * (2π)
        lines!(ax_aniso, kx_b ./ (2π), kz_b ./ (2π); color=:black, linestyle=:dash, linewidth=2)
        text!(ax_aniso, -0.95*Klim_aniso, 0.9*Klim_aniso;
              text=@sprintf("dashed: kᵦ(θ) = kᴿ √|cos θ|,  kᴿ ≈ %.2f", kR),
              color=:black, fontsize=16)
    end

    Colorbar(fig[3, 3], hm_aniso; label="enhancement (decades)")

    Label(fig[0, 1:3], @sprintf("decaying turbulence at t = %.0f", t_m);
          fontsize=26, halign=:center)

    rowgap!(fig.layout, 8)
    colgap!(fig.layout, 12)

    save(outname, fig, px_per_unit=2)
    @info "wrote $outname"
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Usage:
    #   plot_xz_field_and_spectrum.jl [time]
    #   plot_xz_field_and_spectrum.jl <prefix_iso> <prefix_med> [time]
    if length(ARGS) >= 3
        plot_field_and_spectrum(ARGS[1], ARGS[2], parse(Float64, ARGS[3]))
    elseif length(ARGS) == 2
        plot_field_and_spectrum(ARGS[1], ARGS[2], 1000.0)
    else
        tval = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1000.0
        plot_field_and_spectrum("decaying_turbulence_256_9_isotropic",
                                "decaying_turbulence_256_9_medium_surface_waves",
                                tval)
    end
end
