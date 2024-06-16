using GLMakie
using Oceananigans
using Printf

N = 256
name_waves     = "decaying_turbulence_$(N)_surface_waves"
name_isotropic = "decaying_turbulence_$(N)_isotropic"
name_rotating  = "decaying_turbulence_$(N)_rotating"

set_theme!(Theme(fontsize=48))
fig = Figure(resolution=(2800, 960))
n = Observable(1001)

names = [name_rotating, name_waves, name_isotropic]
component = ["ζ", "η", "η"]

for m = 1:3
    name = names[m]
    filename_xy = name * "_xy.jld2"
    filename_xz = name * "_xz.jld2"
    filename_yz = name * "_yz.jld2"

    # wxyt = FieldTimeSeries(filename_xy, "w")
    # wyzt = FieldTimeSeries(filename_yz, "w")
    # wxzt = FieldTimeSeries(filename_xz, "w")

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

    wxy = @lift interior(ωxyt[$n], :, :, 1)
    wyz = @lift interior(ωyzt[$n], 1, :, :)
    wxz = @lift interior(ωxzt[$n], :, 1, :)

    # wxy = @lift interior(wxyt[$n], :, :, 1)
    # wyz = @lift interior(wyzt[$n], 1, :, :)
    # wxz = @lift interior(wxzt[$n], :, 1, :)

    if m == 1
        elevation = 0.7
    else
        elevation = 0.4
    end

    azimuth = 4.0
    perspectiveness = 0.5

    axw = Axis3(fig[2, m]; azimuth, elevation, perspectiveness)
    hidedecorations!(axw)
    hidespines!(axw)

    colormap = :balance
    colorrange = (-2e-2, 2e-2)
    surface!(axw, x_xz, y_xz, z_xz, color=wxz; colormap, colorrange)
    surface!(axw, x_yz, y_yz, z_yz, color=wyz; colormap, colorrange)
    global sf = surface!(axw, x_xy, y_xy, z_xy, color=wxy; colormap, colorrange)
end

title = @lift @sprintf("Vorticity at t = %d", t[$n])
Colorbar(fig[1, 1:3], sf, label=title, vertical=false, flipaxis=true, width=Relative(0.8), height=30)

rowgap!(fig.layout, 1, Relative(-0.1))

display(fig)

record(fig, "decaying_turbulence_$(N)_3d.mp4", 1:Nt, framerate=24) do nn
    @info "Plotting frame $nn of $Nt..."
    n[] = nn
end

# save("decaying_turbulence_3d.png", fig, px_per_unit=5)

