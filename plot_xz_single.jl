#####
##### Single-case xz field-and-spectrum plot: real-space vorticity (η) on top,
##### (kx, kz) energy spectrum on bottom. No isotropic comparison panel.
#####
##### Usage:
#####   julia --project=runs plot_xz_single.jl <prefix> [time]
#####
##### The Stokes shear `∂z_uˢ` for kᵦ overlay is parsed from the prefix: any of
##### {weak, medium, strong, very_strong, very_weak} _surface_waves are supported.
#####

include("plot_xz_spectrum_compare.jl")

using CairoMakie

# ∂z uˢ → ∂²z uˢ for shallow-water Stokes profiles
const STOKES_SHEAR = Dict(
    "very_weak_surface_waves"   => 0.1,
    "weak_surface_waves"        => 0.25,
    "medium_surface_waves"      => 0.5,
    "strong_surface_waves"      => 1.0,
    "very_strong_surface_waves" => 2.0,
)

shear_for_prefix(prefix) = let m = match(r"(very_)?(weak|medium|strong)_surface_waves", prefix)
    m === nothing ? 0.0 : STOKES_SHEAR[m.match]
end

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

function plot_single(prefix, time_target;
                     outname=@sprintf("%s_field_and_spectrum_t%05d.png", prefix, round(Int, time_target)),
                     Klim_aniso=10, log_decades=4, ηpercentile=0.99)

    η, t = eta_xz_at(prefix, time_target)
    Nx, Nz = size(η)
    x = range(0, 1, length=Nx); z = range(0, 1, length=Nz)
    ηlim = quantile(vec(abs.(η)), ηpercentile)
    @info "η colorrange = ±$(round(ηlim, sigdigits=3))"

    kx, kz, E, _, source = spectrum_from_3d_field(prefix, time_target)
    @info "spectrum source: $source"

    vmax = log10(maximum(E)); vmin = vmax - log_decades
    @info "spectrum log range: [$(round(vmin, digits=2)), $(round(vmax, digits=2))]"

    set_theme!(Theme(fontsize=22))
    fig = Figure(size=(1300, 1300))

    # Real-space η
    ax_r = Axis(fig[1, 1]; aspect=1, title=@sprintf("η(x, z)   t = %.0f", t),
                xlabel="x", ylabel="z", xticks=[0, 0.5, 1], yticks=[0, 0.5, 1])
    hm_r = heatmap!(ax_r, x, z, η; colormap=:balance, colorrange=(-ηlim, ηlim))
    Colorbar(fig[1, 2], hm_r; label="η = ∂ₓw − ∂_z u")

    # 2D spectrum (zoomed to barbell scale)
    ax_s = Axis(fig[2, 1]; aspect=1, title="log₁₀ E(kx, kz)",
                xlabel="kx / 2π", ylabel="kz / 2π")
    hm_s = heatmap!(ax_s, kx ./ (2π), kz ./ (2π), log10.(max.(E, 1e-30));
                    colormap=:viridis, colorrange=(vmin, vmax), interpolate=true)
    xlims!(ax_s, -Klim_aniso, Klim_aniso); ylims!(ax_s, -Klim_aniso, Klim_aniso)

    # kᵦ peanut overlay
    ∂²uˢ = shear_for_prefix(prefix)
    if ∂²uˢ > 0
        U = sqrt(2 * sum(E))
        if U > 0
            kR = sqrt(∂²uˢ / U) / (2π)
            θ = range(0, 2π, length=400)
            kβ = @. kR * sqrt(abs(cos(θ)))
            lines!(ax_s, kβ .* cos.(θ), kβ .* sin.(θ);
                   color=:red, linestyle=:dash, linewidth=2)
            text!(ax_s, -0.95*Klim_aniso, 0.9*Klim_aniso;
                  text=@sprintf("kᵦ(θ) = kᴿ √|cos θ|, kᴿ ≈ %.2f", kR),
                  color=:red, fontsize=18)
        end
    end

    Colorbar(fig[2, 2], hm_s; label="log₁₀ E(kx, kz)")

    Label(fig[0, 1:2], prefix; fontsize=20, halign=:center)

    rowgap!(fig.layout, 8); colgap!(fig.layout, 12)
    save(outname, fig, px_per_unit=2)
    @info "wrote $outname"
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        error("usage: julia plot_xz_single.jl <prefix> [time]")
    end
    prefix = ARGS[1]
    tval   = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1000.0
    plot_single(prefix, tval)
end
