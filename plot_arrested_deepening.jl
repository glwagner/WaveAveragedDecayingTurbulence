using Statistics

#=
using GLMakie
using Oceananigans
using Printf

name = "arrested_deepening_wind_strong_waves"
filename_xy = "$(name)_xy.jld2"
filename_yz = "$(name)_yz.jld2"
filename_xz = "$(name)_xz.jld2"
filename_statistics = "$(name)_statistics.jld2"

wwxyt = FieldTimeSeries(filename_xy, "w")
wwyzt = FieldTimeSeries(filename_yz, "w")
wwxzt = FieldTimeSeries(filename_xz, "w")
Bwt = FieldTimeSeries(filename_statistics, "b")
Uwt = FieldTimeSeries(filename_statistics, "u")
W²wt = FieldTimeSeries(filename_statistics, "w²")

name = "arrested_deepening_wind_only"
filename_xy = "$(name)_xy.jld2"
filename_yz = "$(name)_yz.jld2"
filename_xz = "$(name)_xz.jld2"
filename_statistics = "$(name)_statistics.jld2"

wnxyt = FieldTimeSeries(filename_xy, "w")
wnyzt = FieldTimeSeries(filename_yz, "w")
wnxzt = FieldTimeSeries(filename_xz, "w")
Bnt = FieldTimeSeries(filename_statistics, "b")
Unt = FieldTimeSeries(filename_statistics, "u")
W²nt = FieldTimeSeries(filename_statistics, "w²")
=#

grid = wwxyt.grid

tw = wwxyt.times
Nwt = length(tw)

tn = wnxyt.times
Nnt = length(tn)
Nt = min(Nwt, Nnt)
d = Nnt - Nwt
δ = tw[end] / hour - 48

Nx, Ny, Nz = size(grid)
x, y, z = nodes(wnxyt)
Lx = grid.Lx
Ly = grid.Ly
Lz = grid.Lz

W²nzt = zeros(Nz-1, Nt)
W²wzt = zeros(Nz-1, Nt)

∫w²n = zeros(Nt)
∫w²w = zeros(Nt)
Δz = Lz / Nz

for w = 1:Nt
    n = w + d
    W²nzt[:, w] .= interior(W²nt[n], 1, 1, 2:Nz) 
    W²wzt[:, w] .= interior(W²wt[w], 1, 1, 2:Nz) 
    ∫w²n[w] = mean(W²nzt[:, w])
    ∫w²w[w] = mean(W²wzt[:, w])
end

x = collect(x)[:]
y = collect(y)[:]
z = collect(z)[:]

Nz = Nz + 1
x_xz = repeat(x, 1, Nz)
z_xz = repeat(reshape(z, 1, Nz), Nx, 1)
y_xz = 0.001 * ones(Nx, Nz)

y_yz = repeat(y, 1, Nz)
z_yz = repeat(reshape(z, 1, Nz), Ny, 1)
x_yz = 0.001 * ones(Ny, Nz)

# Slight displacements to "stitch" the cube together
x_xy = x
y_xy = y
z_xy = - 0.001 * ones(Nx, Ny)

set_theme!(Theme(fontsize=24, linewidth=3))

fig = Figure(resolution=(1400, 600))

ax1 = Axis3(fig[2, 1]; azimuth=10.7, elevation=0.3, perspectiveness)
ax2 = Axis3(fig[2, 2]; azimuth=10.7, elevation=1.0, perspectiveness)
axU = Axis(fig[2, 3], xlabel="A⁻¹ ∫ u dx dy \n (m s⁻¹)", ylabel="z (m)", xticks=[0, 0.15])
axW = Axis(fig[2, 4], xlabel="A⁻¹ ∫ w² dx dy \n (m² s⁻²)", ylabel="z (m)", yaxisposition=:right, xticks=[0, 0.0004])

hidespines!(axU, :t, :r)
hidespines!(axW, :t, :l)

colgap!(fig.layout, 1, Relative(-0.1))
colgap!(fig.layout, 2, Relative(-0.05))

colsize!(fig.layout, 1, Relative(0.4))
colsize!(fig.layout, 2, Relative(0.4))

colsize!(fig.layout, 3, 150)
colsize!(fig.layout, 4, 150)

#slider = Slider(fig[3, 1:4], range=1:Nwt, startvalue=1)
#w = slider.value
w = Observable(1)
n = @lift $w + d

wwxy = @lift interior(wwxyt[$w], :, :, 1)
wwyz = @lift interior(wwyzt[$w], 1, :, :)
wwxz = @lift interior(wwxzt[$w], :, 1, :)
Uw  = @lift interior(Uwt[$w], 1, 1, :)
W²w = @lift interior(W²wt[$w], 1, 1, :)

wnxy = @lift interior(wnxyt[$n], :, :, 1)
wnyz = @lift interior(wnyzt[$n], 1, :, :)
wnxz = @lift interior(wnxzt[$n], :, 1, :)
Un  = @lift interior(Unt[$n], 1, 1, :)
W²n = @lift interior(W²nt[$n], 1, 1, :)

elevation = 0.4
azimuth = 10.0
perspectiveness = 0.5

title = @lift begin
    t′ = tw[$w] - δ * hour
    if t′ < 0
        @sprintf("t = - %s", prettytime(-t′))
    else
        @sprintf("t = %s", prettytime(t′))
    end
end

Label(fig[1, 1:2], title, tellwidth=false)

hidedecorations!(ax1)
hidedecorations!(ax2)
hidespines!(ax1)
hidespines!(ax2)

colorrange = (-1e-2, 1e-2)
colormap = :balance

surface!(ax1, x_xz, y_xz, z_xz, color=wwxz; colormap, colorrange)
surface!(ax1, x_yz, y_yz, z_yz, color=wwyz; colormap, colorrange)
surface!(ax1, x_xy, y_xy, z_xy, color=wwxy; colormap, colorrange)

surface!(ax2, x_xz, y_xz, z_xz, color=wwxz; colormap, colorrange)
surface!(ax2, x_yz, y_yz, z_yz, color=wwyz; colormap, colorrange)
surface!(ax2, x_xy, y_xy, z_xy, color=wwxy; colormap, colorrange)

grid = Uwt.grid
zc = znodes(grid, Center())
zf = znodes(grid, Face())
lines!(axU, Uw, zc, label="Strong waves")
lines!(axU, Un, zc, label="No waves")
axislegend(axU, position=:rb, labelsize=16)

lines!(axW, W²w, zf, label="Strong waves")
lines!(axW, W²n, zf, label="No waves")

xlims!(axU, -0.01, 0.2)
xlims!(axW, -1e-5, 5e-4)

colsize!(fig.layout, 3, Relative(0.2))
colsize!(fig.layout, 4, Relative(0.2))

display(fig)

record(fig, "arrested_deepening.mp4", 1:Nt, framerate=12) do nn
    @info "Plotting frame $nn of $Nt..."
    w[] = nn
end

using Oceananigans.Units

fig = Figure(resolution=(1140, 840))

axt = Axis(fig[1, 1], yscale=log10, xlabel="Time (hours)", ylabel="⟨w²⟩ (m² s⁻²)", xaxisposition=:top)
lines!(axt, tw ./ hour .- δ, ∫w²n)
lines!(axt, tw ./ hour .- δ, ∫w²w)
ylims!(axt, 5e-8, 3e-5)
xlims!(axt, -3, 24)

axw = Axis(fig[2, 1], xlabel="Time (hours)", ylabel="z (m)")
axn = Axis(fig[3, 1], xlabel="Time (hours)", ylabel="z (m)")

colorrange = (0, 1.5e-4)

heatmap!(axw, tw ./ hour .- δ, zf, W²wzt'; colorrange)
ylims!(axw, -64, 0)
xlims!(axw, -3, 24)

heatmap!(axn, tw ./ hour .- δ, zf, W²nzt'; colorrange)
ylims!(axn, -64, 0)
xlims!(axn, -3, 24)

hidexdecorations!(axw)

display(fig)
