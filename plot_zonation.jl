using GLMakie
using Oceananigans
using Printf

N = 384
medium_waves_name = "decaying_turbulence_$(N)_9_medium_surface_waves"
deep_waves_name = "decaying_turbulence_$(N)_9_deep_surface_waves"
very_deep_waves_name = "decaying_turbulence_$(N)_9_deep_surface_waves"
strong_waves_name = "decaying_turbulence_$(N)_9_strong_surface_waves"

medium_waves_slice_filename    = medium_waves_name * "_xz.jld2"
medium_waves_averages_filename = medium_waves_name * "_averages.jld2"

deep_waves_slice_filename    = deep_waves_name * "_xz.jld2"
deep_waves_averages_filename = deep_waves_name * "_averages.jld2"

very_deep_waves_slice_filename    = very_deep_waves_name * "_xz.jld2"
very_deep_waves_averages_filename = very_deep_waves_name * "_averages.jld2"

strong_waves_slice_filename    = strong_waves_name * "_xz.jld2"
strong_waves_averages_filename = strong_waves_name * "_averages.jld2"

ωmt = FieldTimeSeries(medium_waves_slice_filename, "η")
Umt = FieldTimeSeries(medium_waves_averages_filename, "u")

ωst = FieldTimeSeries(strong_waves_slice_filename, "η")
Ust = FieldTimeSeries(strong_waves_averages_filename, "u")

ωdt = FieldTimeSeries(deep_waves_slice_filename, "η")
Udt = FieldTimeSeries(deep_waves_averages_filename, "u")

ωvt = FieldTimeSeries(very_deep_waves_slice_filename, "η")
Uvt = FieldTimeSeries(very_deep_waves_averages_filename, "u")

t = ωmt.times
Nt = length(t)

set_theme!(Theme(fontsize=24))
fig = Figure(size=(1600, 600))

axo1 = Axis(fig[1, 1], aspect=1, xlabel="x", ylabel="z")
axo2 = Axis(fig[1, 2], aspect=1, xlabel="x", ylabel="z")
axu = Axis(fig[1, 3], xlabel="U", ylabel="z")
axd = Axis(fig[1, 4], xlabel="s", ylabel="z", yaxisposition=:right)

hidespines!(axu, :t, :r)
hidespines!(axd, :l, :t)

hidespines!(axo1)
hidespines!(axo2)

colsize!(fig.layout, 3, Relative(0.3))
colsize!(fig.layout, 4, Relative(0.3))

slider = Slider(fig[2, 1:4], range=1:Nt, startvalue=1)
n = slider.value
ωm = @lift interior(ωmt[$n], :, 1, :)
ωs = @lift interior(ωst[$n], :, 1, :)
ωd = @lift interior(ωdt[$n], :, 1, :)

Um = @lift interior(Umt[$n], 1, 1, :)
Us = @lift interior(Ust[$n], 1, 1, :)
Ud = @lift interior(Udt[$n], 1, 1, :)

levels = -4e-2:1e-3:4e-2
extendhigh = :auto
extendlow = :auto
colorrange = (-5e-2, 5e-2)
colormap = :balance
#kw = (; levels, extendhigh, extendlow, colorrange, colormap)
kw = (; colorrange, colormap)
x, y, z = nodes(ωmt)
z = z .+ 1
#contourf!(axo1, x, z, ωm; kw...)
#contourf!(axo2, x, z, ωs; kw...)
heatmap!(axo1, x, z, ωm; kw...)
heatmap!(axo2, x, z, ωs; kw...)
#Colorbar(fig[1, 1:3], hm; vertical=false, flipaxis=true, label="Vorticity")

z = znodes(Umt) .+ 1
lines!(axu, Um, z)
lines!(axu, Us, z)
lines!(axu, Ud, z)

usm = 0.1 .* z
uss = 0.2 .* z
usd = @. 0.4 * sinh(4 .* z) / sinh(4)
usv = @. 0.8 * exp(8 * (z - 1))

lines!(axd, usm, z)
lines!(axd, uss, z)
lines!(axd, usd, z)

display(fig)

# save("decaying_turbulence.pdf", fig)

