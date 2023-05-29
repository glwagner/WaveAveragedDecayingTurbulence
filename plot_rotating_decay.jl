using GLMakie
using Oceananigans

fig = Figure(resolution=(1200, 600))
ζt = FieldTimeSeries("rotating_decay_slice.jld2", "ζ")
δt = FieldTimeSeries("rotating_decay_slice.jld2", "δ")
t = ζt.times
Nt = length(t)

axz = Axis(fig[1, 1])
axd = Axis(fig[1, 2])
slider = Slider(fig[2, 1:2], range=1:Nt, startvalue=1)
n = slider.value

ζn = @lift interior(ζt[$n], :, :, 1)
δn = @lift interior(δt[$n], :, :, 1)

heatmap!(axz, ζn)
heatmap!(axd, δn)

display(fig)

#=
fig = Figure()
axω = Axis(fig[1, 1], yscale=log10)
axe = Axis(fig[2, 1], yscale=log10)

ξ²t = FieldTimeSeries("rotating_decay_statistics.jld2", "ξ²")
η²t = FieldTimeSeries("rotating_decay_statistics.jld2", "η²")
ζ²t = FieldTimeSeries("rotating_decay_statistics.jld2", "ζ²")

u²t = FieldTimeSeries("rotating_decay_statistics.jld2", "u²")
v²t = FieldTimeSeries("rotating_decay_statistics.jld2", "v²")
w²t = FieldTimeSeries("rotating_decay_statistics.jld2", "w²")
et = FieldTimeSeries("rotating_decay_statistics.jld2", "e")

t = et.times

lines!(axω, t, sqrt.(ξ²t[1, 1, 1, :]), linewidth=4, label="x")
lines!(axω, t, sqrt.(η²t[1, 1, 1, :]), linewidth=4, label="y")
lines!(axω, t, sqrt.(ζ²t[1, 1, 1, :]), linewidth=4, label="z")
axislegend(axω)

lines!(axe, t, u²t[1, 1, 1, :] .* 3/2, label="3u²/2", linewidth=4)
lines!(axe, t, v²t[1, 1, 1, :] .* 3/2, label="3v²/2", linewidth=4)
lines!(axe, t, w²t[1, 1, 1, :] .* 3/2, label="3w²/2", linewidth=4)
lines!(axe, t, et[1, 1, 1, :], label="e", linewidth=4)
axislegend(axe)

display(fig)

=#
