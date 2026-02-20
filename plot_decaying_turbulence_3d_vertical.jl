using GLMakie
using Oceananigans
using Printf

name_waves     = "decaying_turbulence_512_9_medium_surface_waves"
name_isotropic = "decaying_turbulence_512_9_isotropic"
name_rotating  = "decaying_turbulence_512_9_rotating"

set_theme!(Theme(fontsize=48))
fig = Figure(size=(960, 2800))
n = Observable(1)

names = [name_rotating, name_waves, name_isotropic]
component = ["ζ", "η", "η"]

for m = 1:3
    name = names[m]
    filename_xy = name * "_xy.jld2"
    filename_xz = name * "_xz.jld2"
    filename_yz = name * "_yz.jld2"

    ω = component[m]
    ωxyt = FieldTimeSeries(filename_xy, ω)
    ωyzt = FieldTimeSeries(filename_yz, ω)
    ωxzt = FieldTimeSeries(filename_xz, ω)

    global t = ωxyt.times
    global Nt = length(t)

    Nx, Ny, Nz = size(ωxyt.grid)
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

    wxy = @lift Array(interior(ωxyt[$n], :, :, 1))
    wyz = @lift Array(interior(ωyzt[$n], 1, :, :))
    wxz = @lift Array(interior(ωxzt[$n], :, 1, :))

    if m == 1
        elevation = 0.7
    else
        elevation = 0.4
    end

    azimuth = 4.0
    perspectiveness = 0.5

    # Vertical layout: each case in its own row
    axw = Axis3(fig[m, 1]; azimuth, elevation, perspectiveness)
    hidedecorations!(axw)
    hidespines!(axw)

    colormap = :balance
    colorrange = (-15, 15)
    surface!(axw, x_xz, y_xz, z_xz, color=wxz; colormap, colorrange)
    surface!(axw, x_yz, y_yz, z_yz, color=wyz; colormap, colorrange)
    global sf = surface!(axw, x_xy, y_xy, z_xy, color=wxy; colormap, colorrange)
end

display(fig)

record(fig, "decaying_turbulence_512_3d_vertical.mp4", 1:Nt, framerate=12) do nn
    @info "Plotting frame $nn of $Nt..."
    n[] = nn
end

# save("decaying_turbulence_512_3d_vertical.png", fig, px_per_unit=5)
