using GLMakie #CairoMakie
using Oceananigans

filename_long_statistics   = "wavy_kato_phillips_long_waves_statistics.jld2"
filename_none_statistics   = "wavy_kato_phillips_none_waves_statistics.jld2"
filename_medium_statistics = "wavy_kato_phillips_medium_waves_statistics.jld2"
filename_strong_statistics = "wavy_kato_phillips_strong_waves_statistics.jld2"

filename_none_xz   = "wavy_kato_phillips_none_waves_xz.jld2"
filename_long_xz   = "wavy_kato_phillips_long_waves_xz.jld2"
filename_medium_xz = "wavy_kato_phillips_medium_waves_xz.jld2"
filename_strong_xz = "wavy_kato_phillips_strong_waves_xz.jld2"

fig = Figure(resolution=(1200, 400))

bnt = FieldTimeSeries(filename_none_statistics, "b")
blt = FieldTimeSeries(filename_long_statistics, "b")
bmt = FieldTimeSeries(filename_medium_statistics, "b")
bst = FieldTimeSeries(filename_strong_statistics, "b")

wxznt = FieldTimeSeries(filename_none_xz, "w")
wxzlt = FieldTimeSeries(filename_long_xz, "w")
wxzmt = FieldTimeSeries(filename_medium_xz, "w")
wxzst = FieldTimeSeries(filename_strong_xz, "w")

t = blt.times
grid = blt.grid
z = znodes(grid, Center())
Nt = length(t)

axwm = Axis(fig[1, 1])
axwl = Axis(fig[1, 2])
axb = Axis(fig[1, 3])

n = Nt
bn = interior(bnt[n], 1, 1, :)
bl = interior(blt[n], 1, 1, :)
bm = interior(bmt[n], 1, 1, :)
bs = interior(bst[n], 1, 1, :)

lines!(axb, bl, z, label="long")
lines!(axb, bm, z, linewidth=3, label="medium")
lines!(axb, bs, z, linestyle=:dash, linewidth=2, label="strong")
lines!(axb, bn, z, linewidth=3, label="No waves")

axislegend(axb, position=:lt)

x, y, z = nodes(wxzmt)
wlim = 5e-2
colorrange = (-wlim, wlim)
colormap = :balance
heatmap!(axwm, x, z, interior(wxzmt[n], :, 1, :); colorrange, colormap)
heatmap!(axwl, x, z, interior(wxznt[n], :, 1, :); colorrange, colormap)

display(fig)

#save("wavy_kato_phillips.png", fig)

