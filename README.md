# WaveAveragedDecayingTurbulence

Scripts for generating and plotting large eddy simulations of freely-decaying turbulence beneath surface waves used in the paper:

> Wagner, G. L. and Constantinou, N. C. (2025). Phenomenology of decaying turbulence beneath surface waves. *J. Fluid Mech.*, **1020**, A51, doi:[10.1017/jfm.2025.10649](https://doi.org/10.1017/jfm.2025.10649).

The wave-averaged Navier–Stokes equations are integrated with [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl). Stokes drift enters via the Craik–Leibovich vortex force; in the unforced, freely-decaying setting studied here the resulting dynamics is homeomorphic to β-plane turbulence (with `∂²_z u^S` playing the role of β, in the *xz* plane).

## Setup

```julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`GLMakie` is the default for interactive 3D visualization; `CairoMakie` is included for headless plotting (clusters, CI). The simulation drivers run on GPU via CUDA.

---

## Simulation drivers

| script | what it does |
|---|---|
| `decaying_turbulence_256.jl` | Canonical 256³ wave-averaged decaying turbulence runs. Iterates over a list of `kinds` — `isotropic`, `rotating`, `weak`/`medium`/`strong`/`very_strong`/`very_weak`/`deep`/`very_deep` `_surface_waves` — each starting from a common spinup IC (rms ω = 1000 → 10) and integrated with WENO-9 + RK3 to `stop_time = 10⁴`. Writes statistics, xy/xz/yz slices at logarithmically-spaced times, plus full 3D `_fields` at `t = 100, 1000, 10000`. |
| `decaying_turbulence_384.jl` | Same as above at 384³ — used for the longer simulations behind Fig. 3 and Fig. 4 in WC25. |
| `arrested_deepening.jl` | Stratified mixed-layer "arrested deepening" companion run. |
| `runs/decaying_turbulence_256_wise.jl` | Trimmed 256³ driver. `kinds` list is two cases (isotropic + medium-waves) by default; `N`, `STOP_TIME`, and a `SANITY` mode are read from environment variables for short test jobs. |
| `runs/decaying_turbulence_xz25.jl` | **2.5D** xz analogue: `Flat` y, all three velocity components carried so the Craik–Leibovich vortex force is fully active. Runs on CPU or GPU via `ARCH=CPU/GPU`. Cheap sanity-check that the wave-arrest barbell anisotropy survives in 2D dynamics. |

### Slurm scripts (`runs/`)

| script | what it does |
|---|---|
| `runs/sanity.batch` | 15-min sanity check on `gpu-dev` with `N=128`, `SANITY=1`. Boots Julia, verifies the simulation steps, gives a per-iteration timing for extrapolation. |
| `runs/production.batch` | Production run on `gpu-prod` (exclusive node, 4 h). Runs `decaying_turbulence_256_wise.jl` end-to-end. |

---

## Plot scripts: paper figures

| script | what it produces |
|---|---|
| `plot_decaying_turbulence.jl` | Real-space xz/yz/xy snapshots of vorticity for the canonical cases (paper Fig. 1 / Fig. 2 ingredients). |
| `plot_decaying_turbulence_3d.jl` | Stitched-cube 3D visualization of vorticity (paper Fig. 1, three side-by-side cubes). |
| `plot_decaying_turbulence_3d_vertical.jl` | Vertically-stacked variant: rotating, wave-averaged, and isotropic cubes one above the other; renders an mp4. |
| `plot_uw_cross_sections.jl` | u and w cross-sections in the xz plane at a sequence of times — shows the development of large-scale coherent structure under medium waves. |
| `plot_wave_averaged_evolution.jl` | Cross-wave momentum `v(x, y=0, z)` slices at early and late times for deep/medium/weak surface waves, with the corresponding `U(z) = ∫ u dx dy` profiles (paper Fig. 3 layout). |
| `plot_zonation.jl` | Time-series of horizontally-averaged `U(z)` showing depth-alternating jet formation; medium-waves vs deep-waves comparison. |
| `plot_kinetic_energy_decay.jl` | Decay curves `k(t)/k(0)` for one resolution. |
| `plot_many_kinetic_energy_decay.jl` | All six decay curves (isotropic, rotating, three Stokes-shear strengths, deep waves) together with the two-equation-model fits — paper Fig. 4. |
| `plot_isotropic_kinetic_energy_decay.jl` | Resolution-convergence check: isotropic decay across multiple grid sizes. |
| `plot_rotating_decay.jl` | Rotating-only decay curves with the model fits. |
| `plot_arrested_convection.jl`, `plot_arrested_deepening.jl` | Plots for the arrested-deepening / stratified mixed-layer runs. |
| `plot_deep_waves.jl` | Deep-water Stokes-drift case visualization. |
| `plot_kato_phillips.jl` | Plots for the wavy Kato–Phillips runs. |

### Animations

| script | what it produces |
|---|---|
| `animate_decaying_turbulence.jl` | Single-case time animation of vorticity from xz/xy/yz slice data. |
| `mixing_time_movie.jl` | Wind- and wave-driven mixing animation. |

### Spectra

| script | what it produces |
|---|---|
| `spectra.jl` | 1D radial kinetic energy spectrum `E(|k|)` from a 3D field (averaged over y, z). |

---

## Plot scripts: WISE 2026 / wave-arrest "barbell" diagnostic

The wave-averaged dynamics produces an anisotropic energy spectrum in `(kx, kz)` analogous to the β-plane Rhines barbell — energy concentrates in the kx ≈ 0 column (depth-alternating jets). These scripts compute and visualize that signature.

| script | what it produces |
|---|---|
| `plot_xz_spectrum.jl` | Primitives + single-case 2D `(kx, kz)` spectrum from xz slice data. Tukey-window in z to suppress bounded-z FFT leakage; **no horizontal-mean detrending** (that would erase the kx = 0 jet modes). |
| `plot_xz_spectrum_compare.jl` | Side-by-side `log₁₀ E(kx, kz)` for isotropic and medium-waves at a chosen time, with 1D marginal spectra `E(kx)` and `E(kz)`. Reads from full 3D `_fields.jld2` when available and y-averages the 2D FFT for a smooth spectrum. Marks the wave-arrest scale `kᴿ = √(∂²_z u^S / U) / 2π`. |
| `plot_xz_field_and_spectrum.jl` | 2×2 figure: real-space y-vorticity `η = ∂_x w − ∂_z u` on the top row, `(kx, kz)` spectrum on the bottom row, isotropic vs medium-waves. The clearest single picture of the wave-arrest barbell alongside the corresponding real-space coherent structures. |
| `runs/test_spectrum_synthetic.jl` | Sanity test of the spectrum recipe on a synthetic field with a known anisotropy and a wall-anchored mean trend. Produces `synthetic_spectrum_check.png` — the "bad recipe" panel shows the kz-axis ringing the windowing fixes. |

### Quick workflow (WISE-style)

```bash
# Run the simulation (one slurm job)
sbatch runs/production.batch

# Once it lands the t=1000 fields snapshot, plot:
julia --project=. plot_xz_field_and_spectrum.jl 1000
# → xz_field_and_spectrum_t01000.png

julia --project=. plot_xz_spectrum_compare.jl 1000
# → xz_spectrum_compare_t01000.png

# Cheap CPU cross-check:
ARCH=CPU N=192 STOP_TIME=3000 KIND=medium_surface_waves \
  julia --project=. runs/decaying_turbulence_xz25.jl
julia --project=. plot_xz_spectrum.jl decaying_turbulence_xz25_192_medium_surface_waves 3000
```
