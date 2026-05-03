#####
##### Animate the real-space xz vorticity field and its (kx, kz) energy spectrum
##### side-by-side over the full duration of a simulation.
#####
##### Reads from a 2D xz-slice file (`<prefix>_xz.jld2`) — no 3D field data needed.
##### Optionally compares two cases (isotropic + medium-waves) in a 2x2 layout.
#####
##### Usage:
#####   julia --project=. animate_xz_field_and_spectrum.jl <prefix> [out.mp4]
#####       single-case animation
#####
#####   julia --project=. animate_xz_field_and_spectrum.jl <prefix_iso> <prefix_med> [out.mp4]
#####       two-case animation (top row real-space, bottom row spectrum)
#####
##### The animation loops over every time in the xz file. Frame rate defaults to 12 fps.
#####

include("plot_xz_spectrum.jl")

using CairoMakie

"""
    spectrum_from_xz(η, u, v, w; α=0.25)

2D (kx, kz) energy spectrum from a single y-slice. Sums |û|² + |v̂|² + |ŵ|².
Faster than the y-averaged 3D version since we have only one slice.
"""
function spectrum_from_xz(u_xz, v_xz, w_xz; α=0.25)
    Nx_min = minimum(size(a, 1) for a in (u_xz, v_xz, w_xz))
    Nz_min = minimum(size(a, 2) for a in (u_xz, v_xz, w_xz))
    u = u_xz[1:Nx_min, 1:Nz_min]
    v = v_xz[1:Nx_min, 1:Nz_min]
    w = w_xz[1:Nx_min, 1:Nz_min]
    kx, kz, Eu = xz_spectrum(u; α)
    _,  _,  Ev = xz_spectrum(v; α)
    _,  _,  Ew = xz_spectrum(w; α)
    return kx, kz, Eu .+ Ev .+ Ew
end

function load_full_xz_timeseries(prefix)
    f = prefix * "_xz.jld2"
    ut = FieldTimeSeries(f, "u")
    vt = FieldTimeSeries(f, "v")
    wt = FieldTimeSeries(f, "w")
    ηt = FieldTimeSeries(f, "η")
    return (; u=ut, v=vt, w=wt, η=ηt, t=ut.times)
end

squeeze_y(a) = dropdims(a; dims=2)

"""
    animate_one(prefix; outname, framerate, ...)

Single-case animation: real-space η on the left, log10 E(kx, kz) on the right.
Color ranges and v-limits are picked from the early frames so the colormap stays
stable as the field decays (otherwise late-time frames would look saturated by the
shrinking dynamic range).
"""
function animate_one(prefix; outname=prefix * "_xz_animation.mp4",
                     framerate=12, Klim=24, ηpercentile=0.99, log_decades=5)
    ts = load_full_xz_timeseries(prefix)
    Nt = length(ts.t)
    @info "$prefix: $Nt frames over t = $(round(ts.t[1], digits=2)) … $(round(ts.t[end], digits=2))"

    # Lock the η colorrange using the early-mid times (when the field is still
    # active). 99th-percentile across a representative subset of frames.
    sample_idx = unique(round.(Int, range(1, Nt, length=20)))
    η_samples = Float64[]
    for n in sample_idx
        η = squeeze_y(Array(interior(ts.η[n])))
        append!(η_samples, vec(abs.(η)))
    end
    ηlim = quantile(η_samples, ηpercentile)
    @info "η colorrange: ±$(round(ηlim, sigdigits=3))"

    # Lock spectrum vmax from the t≈100-1000 range too (energy peaks there)
    mid_idx = sample_idx[max(1, length(sample_idx) ÷ 3) : end]
    vmaxes = Float64[]
    for n in mid_idx
        u_xz = squeeze_y(Array(interior(ts.u[n])))
        v_xz = squeeze_y(Array(interior(ts.v[n])))
        w_xz = squeeze_y(Array(interior(ts.w[n])))
        _, _, E = spectrum_from_xz(u_xz, v_xz, w_xz)
        push!(vmaxes, log10(maximum(E)))
    end
    vmax = maximum(vmaxes); vmin = vmax - log_decades
    @info "spectrum log range: [$vmin, $vmax]"

    # Set up axes
    set_theme!(Theme(fontsize=22))
    fig = Figure(size=(1300, 600))
    n = Observable(1)

    η0 = squeeze_y(Array(interior(ts.η[1])))
    Nx, Nz = size(η0)
    x = range(0, 1, length=Nx); z = range(0, 1, length=Nz)

    u0 = squeeze_y(Array(interior(ts.u[1])))
    v0 = squeeze_y(Array(interior(ts.v[1])))
    w0 = squeeze_y(Array(interior(ts.w[1])))
    kx, kz, E0 = spectrum_from_xz(u0, v0, w0)

    η_obs = Observable(η0)
    E_obs = Observable(log10.(max.(E0, 1e-30)))

    ax_r = Axis(fig[1, 1]; aspect=1, title="η = ∂ₓw − ∂_z u",
                xlabel="x", ylabel="z", xticks=[0, 0.5, 1], yticks=[0, 0.5, 1])
    hm_r = heatmap!(ax_r, x, z, η_obs; colormap=:balance, colorrange=(-ηlim, ηlim))
    Colorbar(fig[1, 2], hm_r)

    ax_s = Axis(fig[1, 3]; aspect=1, title="log₁₀ E(kx, kz)",
                xlabel="kx / 2π", ylabel="kz / 2π")
    hm_s = heatmap!(ax_s, kx ./ (2π), kz ./ (2π), E_obs;
                    colormap=:viridis, colorrange=(vmin, vmax))
    xlims!(ax_s, -Klim, Klim); ylims!(ax_s, -Klim, Klim)
    Colorbar(fig[1, 4], hm_s)

    title_label = Label(fig[0, 1:4], @sprintf("t = %.1f", ts.t[1]); fontsize=24)

    record(fig, outname, 1:Nt; framerate) do nn
        η_obs[] = squeeze_y(Array(interior(ts.η[nn])))
        u_xz = squeeze_y(Array(interior(ts.u[nn])))
        v_xz = squeeze_y(Array(interior(ts.v[nn])))
        w_xz = squeeze_y(Array(interior(ts.w[nn])))
        _, _, E = spectrum_from_xz(u_xz, v_xz, w_xz)
        E_obs[] = log10.(max.(E, 1e-30))
        title_label.text[] = @sprintf("t = %.1f   (frame %d / %d)", ts.t[nn], nn, Nt)
    end
    @info "wrote $outname"
    return outname
end

"""
    animate_compare(prefix_iso, prefix_med; outname, ...)

2x2 animation: real-space η on top (iso | medium), spectrum on bottom (iso | medium).
"""
function animate_compare(prefix_iso, prefix_med;
                         outname="xz_compare_animation.mp4",
                         framerate=12, Klim=24, ηpercentile=0.99, log_decades=5)
    ts_i = load_full_xz_timeseries(prefix_iso)
    ts_m = load_full_xz_timeseries(prefix_med)

    Nt = min(length(ts_i.t), length(ts_m.t))
    @info "compare: $Nt frames, t = $(round(ts_i.t[1], digits=2)) … $(round(ts_i.t[Nt], digits=2))"

    sample_idx = unique(round.(Int, range(1, Nt, length=20)))
    η_samples = Float64[]
    for n in sample_idx
        for ts in (ts_i, ts_m)
            η = squeeze_y(Array(interior(ts.η[n])))
            append!(η_samples, vec(abs.(η)))
        end
    end
    ηlim = quantile(η_samples, ηpercentile)

    vmaxes = Float64[]
    for n in sample_idx
        for ts in (ts_i, ts_m)
            u = squeeze_y(Array(interior(ts.u[n])))
            v = squeeze_y(Array(interior(ts.v[n])))
            w = squeeze_y(Array(interior(ts.w[n])))
            _, _, E = spectrum_from_xz(u, v, w)
            push!(vmaxes, log10(maximum(E)))
        end
    end
    vmax = maximum(vmaxes); vmin = vmax - log_decades
    @info "η colorrange: ±$(round(ηlim, sigdigits=3)),  spectrum log: [$vmin, $vmax]"

    set_theme!(Theme(fontsize=22))
    fig = Figure(size=(1500, 1300))

    η0_i = squeeze_y(Array(interior(ts_i.η[1])))
    η0_m = squeeze_y(Array(interior(ts_m.η[1])))
    Nx = min(size(η0_i, 1), size(η0_m, 1))
    Nz = min(size(η0_i, 2), size(η0_m, 2))
    x = range(0, 1, length=Nx); z = range(0, 1, length=Nz)

    function spec_n(ts, n)
        u = squeeze_y(Array(interior(ts.u[n])))
        v = squeeze_y(Array(interior(ts.v[n])))
        w = squeeze_y(Array(interior(ts.w[n])))
        spectrum_from_xz(u, v, w)
    end

    kx, kz, E0_i = spec_n(ts_i, 1)
    _,  _,  E0_m = spec_n(ts_m, 1)

    η_obs_i = Observable(η0_i[1:Nx, 1:Nz])
    η_obs_m = Observable(η0_m[1:Nx, 1:Nz])
    E_obs_i = Observable(log10.(max.(E0_i, 1e-30)))
    E_obs_m = Observable(log10.(max.(E0_m, 1e-30)))

    ax_iso_r = Axis(fig[1, 1]; aspect=1, title="isotropic", ylabel="z",
                    xticks=[0, 0.5, 1], yticks=[0, 0.5, 1])
    hm_r = heatmap!(ax_iso_r, x, z, η_obs_i; colormap=:balance, colorrange=(-ηlim, ηlim))
    hidexdecorations!(ax_iso_r, ticks=false)

    ax_med_r = Axis(fig[1, 2]; aspect=1, title="medium surface waves",
                    xticks=[0, 0.5, 1], yticks=[0, 0.5, 1])
    heatmap!(ax_med_r, x, z, η_obs_m; colormap=:balance, colorrange=(-ηlim, ηlim))
    hidexdecorations!(ax_med_r, ticks=false); hideydecorations!(ax_med_r, ticks=false)

    Colorbar(fig[1, 3], hm_r; label="η")

    ax_iso_s = Axis(fig[2, 1]; aspect=1, ylabel="kz / 2π", xlabel="kx / 2π")
    hm_s = heatmap!(ax_iso_s, kx ./ (2π), kz ./ (2π), E_obs_i;
                    colormap=:viridis, colorrange=(vmin, vmax))
    xlims!(ax_iso_s, -Klim, Klim); ylims!(ax_iso_s, -Klim, Klim)

    ax_med_s = Axis(fig[2, 2]; aspect=1, xlabel="kx / 2π")
    heatmap!(ax_med_s, kx ./ (2π), kz ./ (2π), E_obs_m;
             colormap=:viridis, colorrange=(vmin, vmax))
    xlims!(ax_med_s, -Klim, Klim); ylims!(ax_med_s, -Klim, Klim)
    hideydecorations!(ax_med_s, ticks=false)

    Colorbar(fig[2, 3], hm_s; label="log₁₀ E(kx, kz)")

    title_label = Label(fig[0, 1:3], @sprintf("t = %.1f", ts_i.t[1]); fontsize=26)

    record(fig, outname, 1:Nt; framerate) do nn
        η_i = squeeze_y(Array(interior(ts_i.η[nn])))
        η_m = squeeze_y(Array(interior(ts_m.η[nn])))
        η_obs_i[] = η_i[1:Nx, 1:Nz]
        η_obs_m[] = η_m[1:Nx, 1:Nz]
        _, _, E_i = spec_n(ts_i, nn); E_obs_i[] = log10.(max.(E_i, 1e-30))
        _, _, E_m = spec_n(ts_m, nn); E_obs_m[] = log10.(max.(E_m, 1e-30))
        title_label.text[] = @sprintf("t = %.1f   (frame %d / %d)", ts_i.t[nn], nn, Nt)
    end
    @info "wrote $outname"
    return outname
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 0
        # default: compare iso + medium with our production prefixes
        animate_compare("decaying_turbulence_256_9_isotropic",
                        "decaying_turbulence_256_9_medium_surface_waves")
    elseif length(ARGS) == 1
        animate_one(ARGS[1])
    elseif length(ARGS) == 2
        # could be (prefix_iso, prefix_med) or (prefix, outname)
        if endswith(ARGS[2], ".mp4")
            animate_one(ARGS[1]; outname=ARGS[2])
        else
            animate_compare(ARGS[1], ARGS[2])
        end
    else
        animate_compare(ARGS[1], ARGS[2]; outname=ARGS[3])
    end
end
