using GLMakie
using Oceananigans
using Printf

filename = "decaying_turbulence_512_9_medium_surface_waves_xz.jld2"

ut = FieldTimeSeries(filename, "u")
wt = FieldTimeSeries(filename, "w")

t = ut.times

Nx, Ny, Nz = size(ut.grid)
x = range(0, 1, length=Nx)
z = range(-1, 0, length=Nz)

# Find indices closest to t = 20, 400, 1000
targets = [20, 200, 1000]
ns = [argmin(abs.(t .- tᵢ)) for tᵢ in targets]

set_theme!(Theme(fontsize=24))
fig = Figure(size=(1200, 700))

colormap = :balance

for (j, n) in enumerate(ns)
    u = interior(ut[n], :, 1, :)
    w = interior(wt[n], :, 1, :)

    ulim = maximum(abs, u)
    wlim = maximum(abs, w)

    tstr = @sprintf("t = %d", round(Int, t[n]))

    axu = Axis(fig[1, j]; aspect=1, title=tstr,
               xticks=[0, 0.5, 1], yticks=[-1, -0.5, 0])
    axw = Axis(fig[2, j]; aspect=1,
               xticks=[0, 0.5, 1], yticks=[-1, -0.5, 0])

    hmu = heatmap!(axu, x, z, u; colormap, colorrange=(-ulim, ulim))
    hmw = heatmap!(axw, x, z, w; colormap, colorrange=(-wlim, wlim))

    # Top row: hide x decorations, show y only on left
    hidexdecorations!(axu, ticks=false)
    if j > 1
        hideydecorations!(axu, ticks=false)
    else
        axu.ylabel = "z"
    end

    # Bottom row: show x on bottom, show y only on left
    if j > 1
        hideydecorations!(axw, ticks=false)
    else
        axw.ylabel = "z"
    end
    axw.xlabel = "x"

    if j == 3
        Colorbar(fig[1, 4], hmu; label="u")
        Colorbar(fig[2, 4], hmw; label="w")
    end
end

colgap!(fig.layout, 10)

display(fig)

save("uw_cross_sections.png", fig)
