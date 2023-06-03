using GLMakie
using Oceananigans
using Printf

#name = "decaying_turbulence_512_very_deep_surface_waves"
name = "decaying_turbulence_128_navid_surface_waves"
filename_xy = "$(name)_xy.jld2"
filename_yz = "$(name)_yz.jld2"
filename_xz = "$(name)_xz.jld2"

ηxyt = FieldTimeSeries(filename_xy, "η")
ηyzt = FieldTimeSeries(filename_yz, "η")
ηxzt = FieldTimeSeries(filename_xz, "η")

uxyt = FieldTimeSeries(filename_xy, "u")
uyzt = FieldTimeSeries(filename_yz, "u")
uxzt = FieldTimeSeries(filename_xz, "u")

wxyt = FieldTimeSeries(filename_xy, "w")
wyzt = FieldTimeSeries(filename_yz, "w")
wxzt = FieldTimeSeries(filename_xz, "w")

t = ηxyt.times
Nt = length(t)

Nx, Ny, Nz = size(ηxyt.grid)
Δ = 1 / Nx
x = y = z = Δ/2:Δ:1

x_xz = repeat(x, 1, Nz)
z_xz = repeat(reshape(z, 1, Nz), Nx, 1)
y_xz = 0.001 * ones(Nx, Nz)

y_yz = repeat(y, 1, Nz)
z_yz = repeat(reshape(z, 1, Nz), Ny, 1)
x_yz = 0.001 * ones(Ny, Nz)

# Slight displacements to "stitch" the cube together
x_xy = x
y_xy = y
z_xy = 0.998 * ones(Nx, Ny)

n = Observable(501)

ηxy = @lift interior(ηxyt[$n], :, :, 1)
ηyz = @lift interior(ηyzt[$n], 1, :, :)
ηxz = @lift interior(ηxzt[$n], :, 1, :)

uxy = @lift interior(uxyt[$n], :, :, 1)
uyz = @lift interior(uyzt[$n], 1, :, :)
uxz = @lift interior(uxzt[$n], :, 1, :)

wxy = @lift interior(wxyt[$n], :, :, 1)
wyz = @lift interior(wyzt[$n], 1, :, :)
wxz = @lift interior(wxzt[$n], :, 1, :)

elevation = 0.4
azimuth = 4.0
perspectiveness = 0.5

set_theme!(Theme(fontsize=48))
fig = Figure(resolution=(2800, 1200))
axη = Axis3(fig[2, 1]; azimuth, elevation, perspectiveness)
axu = Axis3(fig[2, 2]; azimuth, elevation, perspectiveness)
axw = Axis3(fig[2, 3]; azimuth, elevation, perspectiveness)

title = @lift @sprintf("t = %d", t[$n])
Label(fig[1, 1:2], title)

hidedecorations!(axη)
hidedecorations!(axu)
hidedecorations!(axw)
hidespines!(axη)
hidespines!(axu)
hidespines!(axw)

colorrange = (-5e-2, 5e-2)
colormap = :balance

surface!(axη, x_xz, y_xz, z_xz, color=ηxz; colormap, colorrange)
surface!(axη, x_yz, y_yz, z_yz, color=ηyz; colormap, colorrange)
surface!(axη, x_xy, y_xy, z_xy, color=ηxy; colormap, colorrange)

colorrange = (-2e-4, 2e-4)
surface!(axu, x_xz, y_xz, z_xz, color=uxz; colormap, colorrange)
surface!(axu, x_yz, y_yz, z_yz, color=uyz; colormap, colorrange)
surface!(axu, x_xy, y_xy, z_xy, color=uxy; colormap, colorrange)

surface!(axw, x_xz, y_xz, z_xz, color=wxz; colormap, colorrange)
surface!(axw, x_yz, y_yz, z_yz, color=wyz; colormap, colorrange)
surface!(axw, x_xy, y_xy, z_xy, color=wxy; colormap, colorrange)

display(fig)

record(fig, "deep_wave_decaying_turbulence.mp4", 1:Nt, framerate=24) do nn
    @info "Plotting frame $nn of $Nt..."
    n[] = nn
end

#save("decaying_turbulence.pdf", fig)

