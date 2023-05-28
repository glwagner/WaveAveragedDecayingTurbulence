using Oceananigans
using GLMakie

filename = "wind_and_wave_mixing_slices.jld2"

wt = FieldTimeSeries(filename, "w")
grid = wt.grid
times = wt.times
Nt = length(times)

fig = Figure()
ax = Axis(fig[2, 1])
slider = Slider(fig[1, 1], range=1:Nt, startvalue=1)
n = slider.value

wn = @lift interior(wt[$n], :, 1, :)
heatmap!(ax, wn, colorrange=(-1e-2, 1e-2), colormap=:balance)

display(fig)

