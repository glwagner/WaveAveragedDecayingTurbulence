#####
##### 2x2 figure: real-space xz vorticity (η = ∂xw - ∂zu) on top,
##### (kx, kz) energy spectrum on bottom, side-by-side for isotropic and medium-waves.
#####
##### Mirrors the layout of WC25 Figure 2 (real-space) plus its spectral counterpart.
#####
##### Usage:
#####   julia --project=. plot_xz_field_and_spectrum.jl [time]
#####

include("plot_xz_spectrum_compare.jl")

using CairoMakie

function plot_field_and_spectrum(prefix_iso, prefix_med, time_target;
                                 outname=@sprintf("xz_field_and_spectrum_t%05d.png", round(Int, time_target)),
                                 Klim=24, ηlim=0.012)

    # Real-space xz slice of η at exactly time_target — derived from full 3D field if
    # available (so the field and spectrum panels show the SAME instant).
    function η_xz_at(prefix, t_target)
        fields_file = prefix * "_fields.jld2"
        if isfile(fields_file)
            ut = FieldTimeSeries(fields_file, "u")
            wt = FieldTimeSeries(fields_file, "w")
            n = argmin(abs.(ut.times .- t_target))
            if abs(ut.times[n] - t_target) < 1.0
                u3 = Array(interior(ut[n])); w3 = Array(interior(wt[n]))
                Ny = size(u3, 2); j = max(1, Ny ÷ 2)
                Nz = min(size(u3, 3), size(w3, 3))
                u_xz = u3[:, j, 1:Nz]; w_xz = w3[:, j, 1:Nz]
                Nx_min = min(size(u_xz, 1), size(w_xz, 1))
                u_xz = u_xz[1:Nx_min, :]; w_xz = w_xz[1:Nx_min, :]
                # finite-difference η on the slice (periodic x, regular z)
                Δx = 1 / Nx_min; Δz = 1 / Nz
                ∂x_w = (circshift(w_xz, (-1, 0)) .- circshift(w_xz, (1, 0))) ./ (2Δx)
                ∂z_u = similar(u_xz); ∂z_u[:, 2:end-1] .= (u_xz[:, 3:end] .- u_xz[:, 1:end-2]) ./ (2Δz)
                ∂z_u[:, 1] .= (u_xz[:, 2] .- u_xz[:, 1]) ./ Δz
                ∂z_u[:, end] .= (u_xz[:, end] .- u_xz[:, end-1]) ./ Δz
                return ∂x_w .- ∂z_u, ut.times[n]
            end
        end
        # Fallback to slice file
        return load_xz_slice(prefix, "η", t_target)
    end

    η_i, t_i = η_xz_at(prefix_iso, time_target)
    η_m, t_m = η_xz_at(prefix_med, time_target)

    Nx_i, Nz_i = size(η_i); Nx_m, Nz_m = size(η_m)
    Nx = min(Nx_i, Nx_m); Nz = min(Nz_i, Nz_m)
    η_i = η_i[1:Nx, 1:Nz]; η_m = η_m[1:Nx, 1:Nz]

    x = range(0, 1, length=Nx); z = range(0, 1, length=Nz)

    # Spectra (y-averaged from full 3D field if available)
    kx_i, kz_i, E_i, _, source_i = spectrum_from_3d_field(prefix_iso, time_target)
    kx_m, kz_m, E_m, _, source_m = spectrum_from_3d_field(prefix_med, time_target)
    @info "isotropic source: $source_i, medium source: $source_m"

    Eall = vcat(vec(E_i), vec(E_m))
    vmax = log10(maximum(Eall)); vmin = vmax - 5

    set_theme!(Theme(fontsize=22))
    fig = Figure(size=(1400, 1200))

    # --- Top row: real-space η(x, z) ---
    ax_iso_r = Axis(fig[1, 1]; aspect=1, title="isotropic",
                    xlabel="x", ylabel="z", xticks=[0, 0.5, 1], yticks=[0, 0.5, 1])
    hm_r = heatmap!(ax_iso_r, x, z, η_i; colormap=:balance, colorrange=(-ηlim, ηlim))

    ax_med_r = Axis(fig[1, 2]; aspect=1, title="medium surface waves",
                    xlabel="x", xticks=[0, 0.5, 1], yticks=[0, 0.5, 1])
    heatmap!(ax_med_r, x, z, η_m; colormap=:balance, colorrange=(-ηlim, ηlim))
    hideydecorations!(ax_med_r, ticks=false)

    Colorbar(fig[1, 3], hm_r; label="η = ∂ₓw − ∂_z u")

    # --- Bottom row: 2D spectrum log10 E(kx, kz) ---
    ax_iso_s = Axis(fig[2, 1]; aspect=1, xlabel="kx / 2π", ylabel="kz / 2π")
    hm_s = heatmap!(ax_iso_s, kx_i ./ (2π), kz_i ./ (2π), log10.(max.(E_i, 1e-30));
                    colormap=:viridis, colorrange=(vmin, vmax))
    xlims!(ax_iso_s, -Klim, Klim); ylims!(ax_iso_s, -Klim, Klim)

    ax_med_s = Axis(fig[2, 2]; aspect=1, xlabel="kx / 2π")
    heatmap!(ax_med_s, kx_m ./ (2π), kz_m ./ (2π), log10.(max.(E_m, 1e-30));
             colormap=:viridis, colorrange=(vmin, vmax))
    xlims!(ax_med_s, -Klim, Klim); ylims!(ax_med_s, -Klim, Klim)
    hideydecorations!(ax_med_s, ticks=false)

    # Wave-arrest scale on medium panel
    ∂²uˢ = 0.5
    U = sqrt(2 * sum(E_m))
    if U > 0
        kR = sqrt(∂²uˢ / U) / (2π)
        lines!(ax_med_s, [-Klim, Klim], [ kR,  kR]; color=:red, linestyle=:dash, linewidth=2)
        lines!(ax_med_s, [-Klim, Klim], [-kR, -kR]; color=:red, linestyle=:dash, linewidth=2)
        text!(ax_med_s, -Klim + 1, kR + 1.0; text=@sprintf("kᴿ ≈ %.1f", kR),
              color=:red, fontsize=18)
    end

    Colorbar(fig[2, 3], hm_s; label="log₁₀ E(kx, kz)")

    # Title with time
    Label(fig[0, 1:3], @sprintf("t = %.0f", t_m); fontsize=26, halign=:center)

    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 12)

    save(outname, fig, px_per_unit=2)
    @info "wrote $outname"
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    tval = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1000.0
    plot_field_and_spectrum("decaying_turbulence_256_9_isotropic",
                            "decaying_turbulence_256_9_medium_surface_waves",
                            tval)
end
