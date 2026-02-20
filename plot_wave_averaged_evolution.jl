#using CairoMakie
using GLMakie
using Oceananigans
using Printf

N = 384

filename_very_weak_waves = "decaying_turbulence_$(N)_9_very_weak_surface_waves_xz.jld2"
filename_deep_waves   = "decaying_turbulence_$(N)_9_deep_surface_waves_xz.jld2"
filename_medium_waves = "decaying_turbulence_$(N)_9_medium_surface_waves_xz.jld2"

averages_filename_very_weak_waves = "decaying_turbulence_$(N)_9_very_weak_surface_waves_averages.jld2"
averages_filename_deep_waves   = "decaying_turbulence_$(N)_9_deep_surface_waves_averages.jld2"
averages_filename_medium_waves = "decaying_turbulence_$(N)_9_medium_surface_waves_averages.jld2"

ωmt = FieldTimeSeries(filename_very_weak_waves, "η")
ωdt = FieldTimeSeries(filename_deep_waves, "η")
ωst = FieldTimeSeries(filename_medium_waves, "η")

vmt = FieldTimeSeries(filename_very_weak_waves, "v")
vdt = FieldTimeSeries(filename_deep_waves, "v")
vst = FieldTimeSeries(filename_medium_waves, "v")

Umt = FieldTimeSeries(averages_filename_very_weak_waves, "u")
Udt = FieldTimeSeries(averages_filename_deep_waves, "u")
Ust = FieldTimeSeries(averages_filename_medium_waves, "u")

t = ωmt.times
Nt = length(t)

@show length(Umt)
@show length(Udt)
@show length(Ust)

n1 = 187
n2 = 290
n3 = 322

set_theme!(Theme(fontsize=24))
fig = Figure(size=(1630, 850))

uticks = [-1e-3, 0, 1e-3]
yticks = [-1.0, -0.75, -0.5, -0.25, 0]
axm1 = Axis(fig[2, 3], aspect=1, xlabel="x", ylabel="y", title="Weak waves")
axs1 = Axis(fig[2, 2], aspect=1, xlabel="x", ylabel="z", title="Medium waves")
axd1 = Axis(fig[2, 1], aspect=1, xlabel="x", ylabel="z", title="Deep waves")
axu1 = Axis(fig[2, 4]; xlabel="∫ u dx dy", ylabel="z", yaxisposition=:right, #xticks=uticks,
            xaxisposition=:top, ytrimspine=true, yticks)

axm2 = Axis(fig[3, 3], aspect=1, xlabel="x")
axs2 = Axis(fig[3, 2], aspect=1, xlabel="x") 
axd2 = Axis(fig[3, 1], aspect=1, xlabel="x")
axu2 = Axis(fig[3, 4]; xlabel="∫ u dx dy", ylabel="z", yaxisposition=:right, #xticks=uticks,
            ytrimspine=true, yticks)

            #=
axm3 = Axis(fig[4, 3], aspect=1, xlabel="x", ylabel="y")
axs3 = Axis(fig[4, 2], aspect=1, xlabel="x", ylabel="z")
axd3 = Axis(fig[4, 1], aspect=1, xlabel="x", ylabel="z")
axu3 = Axis(fig[4, 4]; xlabel="∫ u dx dy", ylabel="z", yaxisposition=:right, #xticks=uticks,
            ytrimspine=true, yticks)
            =#

for ax in (axm1, axs1, axd1) #, axm2, axs2, axd2)
    hidedecorations!(ax)
    hidespines!(ax)
end

for ax in (axm2, axs2, axd2)
    hidedecorations!(ax, label=false)
    hidespines!(ax)
end

#=
for ax in (axm3, axs3, axd3)
    hidespines!(ax)
    hideydecorations!(ax)
    hidexdecorations!(ax, label=false)
end
=#

for ax in (axu1, axu2) #, axu3)
    hidespines!(ax, :t, :l, :b)
    xlims!(ax, -1e-2, 2.5e-2)
    #hidexdecorations!(ax, label=false)
end

hidexdecorations!(axu1, grid=false)
hidexdecorations!(axu2, grid=false, label=false)

levels = -4e-2:1e-3:4e-2
extendhigh = :auto
extendlow = :auto
#colorrange = (-3e-2, 3e-2)
colorrange = (-0.01, 0.01)
colormap = :balance
uϵ = 0.01
#kw = (; levels, extendhigh, extendlow, colorrange, colormap)
kw = (; colorrange, colormap)
plotter = heatmap!

lineskw = (; linewidth=6)
colors = Makie.wong_colors(0.6)

z = znodes(Umt)
z = collect(z)

n = n1 
vm = interior(vmt[n], :, 1, :)
vd = interior(vdt[n], :, 1, :)
vs = interior(vst[n], :, 1, :)
heatmap!(axm1, vm; kw...)
heatmap!(axs1, vs; kw...)
hm = heatmap!(axd1, vd; kw...)

Um = interior(Umt[n], 1, 1, :)
Ud = interior(Udt[n], 1, 1, :)
Us = interior(Ust[n], 1, 1, :)
lines!(axu1, Ud .+ 0uϵ, z; color=colors[3], label=" Deep \n waves", lineskw...)
lines!(axu1, Us .+ 1uϵ, z; color=colors[2], label=" Strong \n waves", lineskw...)
lines!(axu1, Um .+ 2uϵ, z; color=colors[1], label=" Medium \n waves", lineskw...)

#vlines!(axu1, 0uϵ, z; color=(:black, 0.5), linewidth=0.5)
#vlines!(axu1, 1uϵ, z; color=(:black, 0.5), linewidth=0.5)
#vlines!(axu1, 2uϵ, z; color=(:black, 0.5), linewidth=0.5)

uˢd = @. uϵ/2 * exp(8 * z)
# lines!(axu1, uˢd, z; color=(:black, 0.5), linewidth=1.0, linestyle=:dash)

n = n2
vm = interior(vmt[n], :, 1, :)
vd = interior(vdt[n], :, 1, :)
vs = interior(vst[n], :, 1, :)
heatmap!(axm2, vm; kw...)
heatmap!(axs2, vs; kw...)
heatmap!(axd2, vd; kw...)

Um = interior(Umt[n], 1, 1, :)
Ud = interior(Udt[n], 1, 1, :)
Us = interior(Ust[n], 1, 1, :)
lines!(axu2, Ud .+ 0uϵ, z; color=colors[3], label=" Deep \n waves", lineskw...)
lines!(axu2, Us .+ 1uϵ, z; color=colors[2], label=" Strong \n waves", lineskw...)
lines!(axu2, Um .+ 2uϵ, z; color=colors[1], label=" Medium \n waves", lineskw...)
#axislegend(axu2, position=:lb, framevisible=false, labelsize=32)

# vlines!(axu2, 0uϵ, z; color=(:black, 0.5), linewidth=0.5)
# vlines!(axu2, 1uϵ, z; color=(:black, 0.5), linewidth=0.5)
# vlines!(axu2, 2uϵ, z; color=(:black, 0.5), linewidth=0.5)

uˢd = @. uϵ/2 * exp(8 * z)
# lines!(axu2, uˢd, z; color=(:black, 0.5), linewidth=1.0, linestyle=:dash)

#=
n = n3
ωm = interior(ωmt[n], :, 1, :)
ωd = interior(ωdt[n], :, 1, :)
ωs = interior(ωst[n], :, 1, :)
heatmap!(axm3, ωm; kw...)
heatmap!(axs3, ωs; kw...)
hm = heatmap!(axd3, ωd; kw...)

Um = interior(Umt[n], 1, 1, :)
Ud = interior(Udt[n], 1, 1, :)
Us = interior(Ust[n], 1, 1, :)
lines!(axu3, Ud .+ 0uϵ, z; color=colors[3], lineskw...)
lines!(axu3, Us .+ 1uϵ, z; color=colors[2], lineskw...)
lines!(axu3, Um .+ 2uϵ, z; color=colors[1], lineskw...)

# vlines!(axu3, 0uϵ, z; color=(:black, 0.5), linewidth=0.5)
# vlines!(axu3, 1uϵ, z; color=(:black, 0.5), linewidth=0.5)
# vlines!(axu3, 2uϵ, z; color=(:black, 0.5), linewidth=0.5)

uˢd = @. uϵ/2 * exp(8 * z)
# lines!(axu3, uˢd, z; color=(:black, 0.5), linewidth=1.0, linestyle=:dash)
=#

Colorbar(fig[1, 1:3], hm; vertical=false, flipaxis=true, label="v(x, y=0, z)")

Label(fig[2, 0], "t = 20", tellheight=false)
Label(fig[3, 0], "t = 400", tellheight=false)
#Label(fig[4, 0], "t = 2×10⁴", tellheight=false)

display(fig)

# save("wave_averaged_evolution.pdf", fig)

