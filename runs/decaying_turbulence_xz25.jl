#####
##### 2.5D x-z analogue of the WC25 wave-averaged decaying turbulence.
#####
##### Domain: x-z (Ny=1, periodic in x, free-slip Bounded in z), but all three velocity
##### components carried, so the Craik-Leibovich vortex force u^S × ω that drives the
##### "barbell" anisotropy in (kx, kz) is fully active. Cheap CPU run for cross-check
##### against the 3D production simulation.
#####
##### Same Stokes shear ∂z u^S = z/2 ("medium_surface_waves") and same initial spectrum
##### shape as the 3D case, scaled to rms ω = 1.
#####

using FFTW
using Oceananigans
using Statistics
using Printf
using JLD2

using Oceananigans.Operators

@inline ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2

@inline function tke25(i, j, k, grid, u, v, w)
    u² = ℑxᶜᵃᵃ(i, j, k, grid, ϕ², u)
    v² = ℑyᵃᶜᵃ(i, j, k, grid, ϕ², v)
    w² = ℑzᵃᵃᶜ(i, j, k, grid, ϕ², w)
    return (u² + v² + w²) / 2
end

@inline y_enstrophy25(i, j, k, grid, u, v, w) = (∂xᶠᶜᶠ(i, j, k, grid, w) - ∂zᶠᶜᶠ(i, j, k, grid, u))^2

const N = parse(Int,     get(ENV, "N",         "256"))
const STOP_TIME = parse(Float64, get(ENV, "STOP_TIME", "1e4"))
const KIND      = get(ENV, "KIND", "medium_surface_waves")  # or "isotropic"

Nx = N; Nz = N
weno_order = 9
cfl = 0.5
arch = get(ENV, "ARCH", "GPU") == "GPU" ? Oceananigans.GPU() : Oceananigans.CPU()

x = (0, 1); z = (0, 1)
topology = (Periodic, Flat, Bounded)
grid = RectilinearGrid(arch, size=(Nx, Nz), halo=(7, 7); x, z, topology)

# 2D random IC with the same kx,kz target spectrum as the 3D run, then scale rms ω to 1.
function velocity_spectral_shape(k)
    k₀ = 32 * 2π
    k′ = k / k₀
    return k * sqrt(exp(-2 * k′^2))
end

function initial_xz_field(Nx, Nz)
    kx = vcat(0:Nx÷2, -Nx÷2+1:-1) .* (2π)
    kz = vcat(0:Nz÷2, -Nz÷2+1:-1) .* (2π)
    KX = repeat(reshape(kx, Nx, 1), 1, Nz)
    KZ = repeat(reshape(kz, 1, Nz), Nx, 1)
    K  = sqrt.(KX.^2 .+ KZ.^2)
    f̂  = velocity_spectral_shape.(K) .* (randn(Nx, Nz) .+ im .* randn(Nx, Nz))
    f̂[1, 1] = 0
    return real.(ifft(f̂))
end

import Random
Random.seed!(20260501)
u0_xz = initial_xz_field(Nx, Nz)
v0_xz = initial_xz_field(Nx, Nz)
w0_xz = initial_xz_field(Nx, Nz)

# Flat y → arrays are (Nx, Nz)
u0 = u0_xz
v0 = v0_xz
w0 = zeros(Nx, Nz+1); w0[:, 2:Nz] .= w0_xz[:, 2:Nz]

struct ShallowStokesShear; shear::Float64; end
@inline (s::ShallowStokesShear)(z, t) = s.shear * z

if KIND == "medium_surface_waves"
    kwargs = (; stokes_drift = UniformStokesDrift(∂z_uˢ=ShallowStokesShear(0.5)))
elseif KIND == "isotropic"
    kwargs = NamedTuple()
else
    error("unknown KIND=$KIND")
end

timestepper = :RungeKutta3
advection = WENO(order=weno_order)

model = NonhydrostaticModel(; grid, timestepper, advection, kwargs...)
set!(model, u=u0, v=v0, w=w0)

# Subtract mean and scale to rms(ω) = 1
u, v, w = model.velocities
parent(u) .-= mean(u); parent(v) .-= mean(v); parent(w) .-= mean(w)

η² = KernelFunctionOperation{Face, Center, Face}(y_enstrophy25, grid, u, v, w)
η² = compute!(Field(η²))
ω0 = sqrt(mean(η²))
parent(u) ./= ω0; parent(v) ./= ω0; parent(w) ./= ω0

max_u = maximum(abs, model.velocities.u)
Δx = 1 / Nx
Δt = 1e-3 * Δx / max_u
simulation = Simulation(model; Δt, stop_time=STOP_TIME)
conjure_time_step_wizard!(simulation, IterationInterval(3); cfl)

prefix = "decaying_turbulence_xz25_$(Nx)_$(KIND)"

η = ∂x(w) - ∂z(u)
e  = KernelFunctionOperation{Center, Center, Center}(tke25, grid, u, v, w)
e  = Field(e); η² = Field(η²)
E  = Field(Average(e)); Y² = Field(Average(η²))

start_time = Ref(time_ns())
function progress(sim)
    msg = @sprintf("Iter: % 6d, t: %7.2f, Δt: %.2e", iteration(sim), time(sim), sim.Δt)
    u, v, w = sim.model.velocities
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e)", maximum(abs, u), maximum(abs, v), maximum(abs, w))
    compute!(e); compute!(η²)
    msg *= @sprintf(", ⟨e⟩: %.2e, rms(η): %.2e", mean(e), sqrt(mean(η²)))
    elapsed = (time_ns() - start_time[]) / 1e9
    msg *= @sprintf(", wall: %s", prettytime(elapsed))
    start_time[] = time_ns()
    @info msg
end
simulation.callbacks[:progress] = Callback(progress, IterationInterval(200))

slices_saves = 200
n_slices = range(-1, stop=log10(STOP_TIME), length=slices_saves)
slices_times = 10 .^ n_slices

simulation.output_writers[:xz] = JLD2OutputWriter(model, (; η, u, v, w, e);
                                                  schedule = SpecifiedTimes(slices_times),
                                                  filename = prefix * "_xz",
                                                  overwrite_existing = true)

simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w);
                                                     schedule = SpecifiedTimes(100.0, 1000.0, Float64(STOP_TIME)),
                                                     filename = prefix * "_fields",
                                                     with_halos = true,
                                                     overwrite_existing = true)

@info "Running 2.5D xz $(KIND), N=$N, stop_time=$STOP_TIME..."
run!(simulation)
