using GLMakie
using Oceananigans
using Printf

N = 256
filename_waves     = "decaying_turbulence_$(N)_surface_waves_slice.jld2"
filename_isotropic = "decaying_turbulence_$(N)_isotropic_slice.jld2"
filename_rotating  = "decaying_turbulence_$(N)_rotating_slice.jld2"

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

n = Observable(1)
ωw = @lift interior(ωwt[$n], :, 1, :)
ωi = @lift interior(ωit[$n], :, 1, :)
ωr = @lift interior(ωrt[$n], :, :, 1)

colorrange = (-5e-2, 5e-2)
hm = heatmap!(axr, ωr; colorrange, colormap=:balance)

label = @lift @sprintf("Vorticity at t = % 5d", t[$n])
Colorbar(fig[1, 1:3], hm; vertical=false, flipaxis=true, label)

heatmap!(axw, ωw; colorrange, colormap=:balance)
heatmap!(axi, ωi; colorrange, colormap=:balance)

hidedecorations!(axi, label=false)
hidedecorations!(axr, label=false)
hidedecorations!(axw, label=false)

display(fig)

record(fig, "decaying_turbulence_$(N)_visualization.mp4", 1:Nt, framerate=24) do nn
    @info "Making frame $nn of $Nt..."
    n[] = nn
end

