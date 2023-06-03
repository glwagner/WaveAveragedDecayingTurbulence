using CairoMakie
using Oceananigans
using Printf

filename_waves     = "decaying_turbulence_256_surface_waves_slice.jld2"
filename_isotropic = "decaying_turbulence_256_isotropic_slice.jld2"
filename_rotating  = "decaying_turbulence_256_rotating_slice.jld2"

set_theme!(Theme(fontsize=48))

fig = Figure(resolution=(1600, 700))

ωwt = FieldTimeSeries(filename_waves, "η")
ωit = FieldTimeSeries(filename_isotropic, "η")
ωrt = FieldTimeSeries(filename_rotating, "ζ")

t = ωwt.times
Nt = length(t)

axr = Axis(fig[2, 1], aspect=1, xlabel="x", ylabel="y")
axw = Axis(fig[2, 2], aspect=1, xlabel="x", ylabel="z")
axi = Axis(fig[2, 3], aspect=1, xlabel="x", ylabel="z")

hidedecorations!(axr)
hidedecorations!(axw)
hidedecorations!(axi)
hidespines!(axr)
hidespines!(axw)
hidespines!(axi)

n = Observable(501)
ωw = @lift interior(ωwt[$n], :, 1, :)
ωi = @lift interior(ωit[$n], :, 1, :)
ωr = @lift interior(ωrt[$n], :, :, 1)

levels = -4e-2:1e-3:4e-2
extendhigh = :auto
extendlow = :auto
colorrange = (-5e-2, 5e-2)
colormap = :balance
kw = (; levels, extendhigh, extendlow, colorrange, colormap)
hm = contourf!(axr, ωr; kw...)

Colorbar(fig[1, 1:3], hm; vertical=false, flipaxis=true, label="Vorticity")

contourf!(axw, ωw; kw...)
contourf!(axi, ωi; kw...)

hidedecorations!(axi, label=false)
hidedecorations!(axr, label=false)
hidedecorations!(axw, label=false)

display(fig)

save("decaying_turbulence.pdf", fig)

