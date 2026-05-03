# WaveAveragedDecayingTurbulence

Scripts for generating and plotting large eddy simulations of freely-decaying turbulence beneath surface waves used in the paper:

> Wagner, G. L. and Constantinou, N. C. (2025). Phenomenology of decaying turbulence beneath surface waves. *J. Fluid Mech.*, **1020**, A51, doi:[10.1017/jfm.2025.10649](https://doi.org/10.1017/jfm.2025.10649). [Open-access PDF.](https://glwagner.github.io/assets/pdf/wave_averaged_decaying_turbulence.pdf)

The wave-averaged Navier–Stokes equations are integrated with [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl). Stokes drift enters via the Craik–Leibovich vortex force; in the unforced, freely-decaying setting studied here the resulting dynamics is **homeomorphic to β-plane turbulence** (with `∂²_z u^S` playing the role of β, in the *xz* plane instead of *xy*). The signature feature is a **wave-arrest "barbell"** in the (kx, kz) energy spectrum — the analog of the Rhines dumbbell — and depth-alternating zonal jets in real space.

---

## Table of contents

- [Background reading](#background-reading)
- [Installing Julia](#installing-julia)
- [Setting up this project](#setting-up-this-project)
- [Quick start](#quick-start)
- [Simulation drivers](#simulation-drivers)
- [Slurm batch scripts](#slurm-batch-scripts)
- [Plot scripts: paper figures](#plot-scripts-paper-figures)
- [Plot scripts: wave-arrest "barbell" diagnostic](#plot-scripts-wave-arrest-barbell-diagnostic)
- [Animations](#animations)
- [Glossary of variables and quantities](#glossary-of-variables-and-quantities)

---

## Background reading

If you're new to the wave-averaged Navier–Stokes equations, β-plane turbulence, or the simulation tooling, these are the canonical references:

**This paper (the one this repo accompanies):**
- Wagner, G. L. & Constantinou, N. C. 2025. *Phenomenology of decaying turbulence beneath surface waves*. JFM **1020**, A51. [Open-access PDF.](https://glwagner.github.io/assets/pdf/wave_averaged_decaying_turbulence.pdf)

**Wave-averaged Navier–Stokes (Craik–Leibovich):**
- Craik, A. D. D. & Leibovich, S. 1976. *A rational model for Langmuir circulations*. JFM **73**, 401–426.
- Holm, D. D. 1996. *The ideal Craik–Leibovich equations*. *Physica D* **98**, 415–441. — derives the vortex force interpretation and the rotating-turbulence analog.
- Suzuki, N. & Fox-Kemper, B. 2016. *Understanding Stokes forces in the wave-averaged equations*. JGR **121**, 3579–3596.
- van den Bremer, T. & Breivik, Ø. 2018. *Stokes drift*. *Phil. Trans. R. Soc. A* **376**, 20170104.

**β-plane turbulence and the Rhines / Vallis–Maltrud barbell:**
- Rhines, P. B. 1975. *Waves and turbulence on a beta-plane*. JFM **69**, 417–443. — the seminal Rhines-arrest paper, defines `kᵦ = √(β/U)`.
- Vallis, G. K. & Maltrud, M. E. 1993. *Generation of mean flows and jets on a beta plane and over topography*. JPO **23**, 1346–1362.
- Vallis, G. K. 2017. *Atmospheric and Oceanic Fluid Dynamics*, 2nd ed., Cambridge UP — Chapter 12, especially Fig 12.3, shows the canonical dumbbell evolution.

**Oceananigans.jl (the LES code):**
- [Documentation](https://clima.github.io/OceananigansDocumentation/stable/) — start with the *Quick start* and *Models → NonhydrostaticModel*.
- [GitHub](https://github.com/CliMA/Oceananigans.jl).
- [Discussions](https://github.com/CliMA/Oceananigans.jl/discussions) — best place to ask questions.
- Ramadhan, A. *et al.* 2020. *Oceananigans.jl: Fast and friendly geophysical fluid dynamics on GPUs*. *J. Open Source Soft.* **5**, 2018.
- Wagner, G. L. *et al.* 2025. *High-level, high-resolution ocean modelling at all scales with Oceananigans*. arXiv:2502.14148.

**Julia language (if you've never used it):**
- [Julia documentation](https://docs.julialang.org/en/v1/).
- [Modern Julia Workflows](https://modernjuliaworkflows.org/) — practical guide to project envs, REPL, debugging.

---

## Installing Julia

Use the official installer [`juliaup`](https://github.com/JuliaLang/juliaup), which manages multiple Julia versions cleanly:

**Linux / macOS:**
```bash
curl -fsSL https://install.julialang.org | sh
```

**Windows:**
```powershell
winget install julia -s msstore
```

After installation, restart your shell and verify:
```bash
julia --version
```

This project is currently tested on **Julia 1.10**. To pin that version specifically:
```bash
juliaup add 1.10
juliaup default 1.10
```

(Newer Julia versions also work; some packages may emit a "tested on 1.10" warning.)

---

## Setting up this project

Clone the repo and instantiate the Julia environment from the committed `Manifest.toml` (this resolves the exact package versions the paper used):

```bash
git clone https://github.com/glwagner/WaveAveragedDecayingTurbulence.git
cd WaveAveragedDecayingTurbulence
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

First-time instantiation downloads/precompiles Oceananigans, CUDA, FFTW, Makie, etc., and may take 10–30 minutes.

**GPU note:** `CUDA.jl` requires an NVIDIA driver compatible with CUDA 11 or 12. If your cluster has a CUDA 13 driver, the published `Manifest.toml` won't work — use the secondary `runs/Project.toml` env instead, which pins Oceananigans `0.95.x` but resolves to a newer CUDA.jl:
```bash
julia --project=runs -e 'using Pkg; Pkg.instantiate()'
```

**Headless plotting:** the project includes both `GLMakie` (interactive, needs an X11/OpenGL display) and `CairoMakie` (file-only, works on a bare cluster). Use `CairoMakie` on headless machines — most plot scripts in this repo can be switched by editing `using GLMakie` → `using CairoMakie`. The newer barbell-diagnostic scripts (`plot_xz_*`, `animate_xz_*`) already use `CairoMakie`.

---

## Quick start

```bash
# 1. Reproduce the paper's 256³ runs (GPU; takes hours; outputs *.jld2)
julia --project=. decaying_turbulence_256.jl

# 2. Reproduce the paper figures
julia --project=. plot_decaying_turbulence_3d.jl   # Fig 1
julia --project=. plot_wave_averaged_evolution.jl  # Fig 3
julia --project=. plot_many_kinetic_energy_decay.jl # Fig 4

# 3. Generate the (kx, kz) "barbell" energy spectrum diagnostic
julia --project=. plot_xz_field_and_spectrum.jl <iso_prefix> <med_prefix> 1000
# or single-case
julia --project=. plot_xz_single.jl decaying_turbulence_256_9_medium_surface_waves 1000

# 4. Cheap 2.5D analogue: a *single* xz simulation that captures the same wave
#    arrest physics, runs in ~10 minutes on one GPU at N=256, t=10⁴
KIND=very_strong_surface_waves N=256 julia --project=runs runs/decaying_turbulence_xz25.jl
julia --project=runs plot_xz_single.jl \
    decaying_turbulence_xz25_256_very_strong_surface_waves 3000
```

---

## Simulation drivers

| script | what it does |
|---|---|
| `decaying_turbulence_256.jl` | **Canonical 256³ wave-averaged decaying turbulence runs.** Iterates over a list of `kinds` — `isotropic`, `rotating`, `weak`/`medium`/`strong`/`very_strong`/`very_weak`/`deep`/`very_deep` `_surface_waves` — each starting from a common spinup IC (rms ω = 1000 → 10) and integrated with WENO-9 + RK3 to `stop_time = 10⁴`. Writes statistics, `xy/xz/yz` slices at logarithmically-spaced times, plus full 3D `_fields` at `t = 100, 1000, 10000`. **GPU only.** |
| `decaying_turbulence_384.jl` | Same as above at **384³** — used for the longer simulations behind Fig 3 and Fig 4 of the paper. |
| `arrested_deepening.jl` | Stratified mixed-layer "arrested deepening" companion run (not in the JFM paper). |
| `runs/decaying_turbulence_256_wise.jl` | **Trimmed 256³ driver** for the WISE-conference figures. `kinds` defaults to `medium_surface_waves` only; `N`, `STOP_TIME`, and a `SANITY` mode are read from environment variables for short test jobs. |
| `runs/decaying_turbulence_xz25.jl` | **2.5D xz analogue** (recommended for cheap exploration). Uses Oceananigans' `Flat` y-topology; all three velocity components are carried so the Craik–Leibovich vortex force is fully active. Configurable via env vars: `N`, `STOP_TIME`, `KIND`, `LX`, `LZ`, `ARCH=CPU/GPU`. At `N=256`, `LX=1`, GPU, runs in ~10 min to t=10⁴. |

### Initial condition recipe

All drivers use the same recipe for the IC (described in §2.1 of the paper):

1. Sample velocity in spectral space with shape `|û(k)|² ∝ k² exp(-2 (k/k₀)²)`, `k₀ = 32·2π`.
2. Spin up with the inviscid, Stokes-drift-free model from rms vorticity 1000 down to rms 10 — this rounds out the IC into a quasi-equilibrium turbulent field.
3. Save the spun-up state to `initial_conditions_<N>.jld2`. This file is reused across all `kinds` for a given resolution so the only difference between runs is the wave forcing.

---

## Slurm batch scripts

Located under `runs/`. Adjust partitions and node-list for your cluster.

| script | what it does |
|---|---|
| `runs/sanity.batch` | 15-minute sanity check on `gpu-dev`: instantiates Julia, runs `N=128` for ~5 simulation time units, reports per-iter timing for extrapolation. |
| `runs/production.batch` | Production run on `gpu-prod` (exclusive node, 4 h budget). Calls `runs/decaying_turbulence_256_wise.jl`. |
| `runs/xz25_long.batch` | 2.5D xz analogue on `gpu-prod`. Reads `KIND`, `N`, `LX`, `LZ`, `STOP_TIME`, `ARCH` from env. Submit several in parallel for a parameter sweep. |

Examples:

```bash
# Sanity check first run on a fresh cluster
sbatch runs/sanity.batch

# Full 3D production
sbatch runs/production.batch

# 2.5D parameter sweep over Stokes shear strength
KIND=weak_surface_waves        sbatch runs/xz25_long.batch
KIND=medium_surface_waves      sbatch runs/xz25_long.batch
KIND=strong_surface_waves      sbatch runs/xz25_long.batch
KIND=very_strong_surface_waves sbatch runs/xz25_long.batch

# Larger box for finer spectral resolution
N=512 LX=2 LZ=2 KIND=very_strong_surface_waves sbatch runs/xz25_long.batch
```

---

## Plot scripts: paper figures

These reproduce or are adapted from the figures in Wagner & Constantinou 2025.

| script | what it produces |
|---|---|
| `plot_decaying_turbulence.jl` | Real-space `xz/yz/xy` snapshots of vorticity for the canonical cases (paper Fig 1 / Fig 2 ingredients). |
| `plot_decaying_turbulence_3d.jl` | Stitched-cube 3D visualization of vorticity (paper Fig 1, three side-by-side cubes). |
| `plot_decaying_turbulence_3d_vertical.jl` | Vertically-stacked variant: rotating, wave-averaged, and isotropic cubes one above the other; renders an mp4. |
| `plot_uw_cross_sections.jl` | u and w cross-sections in the xz plane at a sequence of times — shows the development of large-scale coherent structure under medium waves. |
| `plot_wave_averaged_evolution.jl` | Cross-wave momentum `v(x, y=0, z)` slices at early and late times for deep / medium / weak surface waves, with the corresponding `U(z) = ∫ u dx dy` profiles (paper **Fig 3** layout). |
| `plot_zonation.jl` | Time-series of horizontally-averaged `U(z)` showing depth-alternating jet formation; medium-waves vs deep-waves comparison. |
| `plot_kinetic_energy_decay.jl` | KE decay curves `k(t) / k(0)` for one resolution. |
| `plot_many_kinetic_energy_decay.jl` | All six decay curves (isotropic, rotating, three Stokes-shear strengths, deep waves) together with the two-equation-model fits — paper **Fig 4**. |
| `plot_isotropic_kinetic_energy_decay.jl` | Resolution-convergence check: isotropic decay across multiple grid sizes. |
| `plot_rotating_decay.jl` | Rotating-only decay curves with the model fits. |
| `plot_arrested_convection.jl`, `plot_arrested_deepening.jl` | Plots for the arrested-deepening / stratified mixed-layer runs. |
| `plot_deep_waves.jl` | Deep-water Stokes-drift case visualization. |
| `plot_kato_phillips.jl` | Plots for the wavy Kato–Phillips runs. |
| `spectra.jl` | 1D radial kinetic energy spectrum `E(|k|)` from a 3D field (averaged over y, z). |

---

## Plot scripts: wave-arrest "barbell" diagnostic

The wave-averaged dynamics produces an anisotropic energy spectrum in `(kx, kz)` analogous to the **β-plane Rhines barbell** (Vallis 2017, Fig 12.3): energy concentrates in the `kx ≈ 0` column (depth-alternating jets), with a peanut-shaped depleted region around it bounded by the wave-arrest scale `kᴿ = √(∂²_z u^S / U) / 2π`.

| script | what it produces |
|---|---|
| `plot_xz_spectrum.jl` | Primitives + single-case 2D `(kx, kz)` spectrum from xz slice data. Tukey-window in z to suppress bounded-z FFT leakage; **no horizontal-mean detrending** (that would erase the kx = 0 jet modes). Reads box size from the grid for correct wavenumber units. |
| `plot_xz_spectrum_compare.jl` | Side-by-side `log₁₀ E(kx, kz)` for two cases at a chosen time, with 1D marginal spectra `E(kx)` and `E(kz)`. Reads from the full 3D `_fields.jld2` when available and y-averages the 2D FFT for a smooth spectrum. |
| `plot_xz_field_and_spectrum.jl` | 3-panel comparison figure: real-space y-vorticity `η = ∂_x w − ∂_z u` on top, raw `(kx, kz)` spectra in the middle, and the wave-induced anisotropy `log₁₀(E_medium / E_isotropic)` on the bottom — that last panel is where the barbell pops. |
| `plot_xz_single.jl` | **Single-case** version: real-space η on top, log₁₀ E(kx, kz) on bottom, with the kᵦ peanut overlay. Use when you don't want an isotropic baseline. |
| `animate_xz_field_and_spectrum.jl` | Animation of η + spectrum evolving over the full simulation. Single-case (one prefix → one mp4) or two-case comparison (two prefixes → 2×2 mp4). Spectra computed directly from 2D xz slice data — works on any saved xz output without needing 3D fields. |
| `runs/test_spectrum_synthetic.jl` | Sanity test of the spectrum recipe on a synthetic field with a known anisotropy and a wall-anchored mean trend. The "bad recipe" (no detrending, no z-window) panel shows the kz-axis ringing the windowing fixes. |

### Recommended workflow for the barbell

```bash
# 1. Run the 2.5D xz simulation (cheap; one GPU, ~10 min for N=256, t=10⁴).
N=256 KIND=very_strong_surface_waves \
  julia --project=runs runs/decaying_turbulence_xz25.jl

# 2. Plot the field + spectrum at a few times. Pre-collapse times (t≈300–3000)
#    show the most informative dumbbell; t=10⁴ is fully arrested onto the box
#    scale and can look "too clean" (just a couple of saturated cells).
for t in 100 300 1000 3000 10000; do
  julia --project=runs plot_xz_single.jl \
    decaying_turbulence_xz25_256_very_strong_surface_waves $t
done

# 3. (Optional) Larger box → finer spectral resolution. With LX=2 we get Δk = π
#    instead of 2π, so the dumbbell spans ~8 cells from origin to kᵦ instead of ~4.
N=512 LX=2 LZ=2 KIND=very_strong_surface_waves \
  julia --project=runs runs/decaying_turbulence_xz25.jl

# 4. Animation of the cascade developing into the dumbbell:
julia --project=runs animate_xz_field_and_spectrum.jl \
  decaying_turbulence_xz25_512_L2_very_strong_surface_waves
```

### Tuning the spectrum visualization

`plot_xz_single.jl` and `plot_xz_field_and_spectrum.jl` accept the following knobs at the top of the function:
- `Klim_aniso` — half-range of the (kx, kz) plot in 2π units (default 8–10). Zoom further in to make low-k structure visible.
- `log_decades` — dynamic range of the log color scale (default 1.5). Tight is better when the dumbbell hollow is only ~0.5–1 decade below the peak.
- `ηpercentile` — percentile-based color saturation for the real-space η panel (default 0.99), so a few extreme vortices don't blow out the colorscale.

---

## Animations

| script | what it produces |
|---|---|
| `animate_decaying_turbulence.jl` | Single-case time animation of vorticity from `xz/xy/yz` slice data. |
| `mixing_time_movie.jl` | Wind- and wave-driven mixing animation. |
| `animate_xz_field_and_spectrum.jl` | (See barbell diagnostics above.) η + spectrum side-by-side over the full simulation; auto-detects single-case vs two-case mode from the number of prefix arguments. |

---

## Glossary of variables and quantities

For readers coming from outside the wave-averaged-turbulence literature, this is what the symbols in the scripts and figures mean.

- **`u, v, w`** — Lagrangian-mean velocity components in `x, y, z`. The Lagrangian mean `u = u_E + u^S` is the sum of the Eulerian-mean velocity and the Stokes drift correction (van den Bremer & Breivik 2018).
- **`u^S(z)`** — Stokes drift profile. Shallow-water case: `u^S(z) ≈ (1+z²/2)/2`, so `∂_z u^S = z/2` (linear in z) and `∂²_z u^S = 1/2` (constant). The drivers parameterize this as a `ShallowStokesShear(shear)` callable returning `∂_z u^S = shear · z`; `shear = 0.1, 0.25, 0.5, 1.0, 2.0` for `very_weak / weak / medium / strong / very_strong`.
- **`Ω`** — background "vorticity-like" vector that appears in the wave-averaged momentum equation as `Ω × u` (eq 2.1 in Wagner & Constantinou 2025).
  - Rotating turbulence: `Ω = f ẑ` (Coriolis).
  - Wave-averaged turbulence: `Ω = -∇ × u^S = -∂_z u^S ŷ`.
- **`η = ∂_x w − ∂_z u`** — y-component of vorticity. The wave-averaged paper plots this in Fig 1(b), Fig 2(b) — coherent structures form perpendicular to the wave-propagation direction (`x`).
- **`ζ = ∂_x v − ∂_y u`** — z-component of vorticity (the "rotating-turbulence analog" of η).
- **`kᴿ = √(∂²_z u^S / U) / 2π`** — wave-arrest wavenumber, the analog of the Rhines wavenumber `kᵦ = √(β/U)`. Energy concentrates in modes with `k < kᴿ`. Real-space Rhines length is `ℓᴿ = 2π / kᴿ`.
- **`Ps = |∇×u| / |∂_z u^S|`** — the "pseudovorticity number" (§4 of WC25), the analog of the Rossby number for surface-wave turbulence.
- **"Depth-alternating jets" / "zonation"** — the analog of zonal jets in β-plane turbulence: alternating bands of `u(z)` that don't vary much in `x`, i.e., spectral content concentrated at `kx = 0`, finite `kz`.
- **"Barbell" / "dumbbell"** — the peanut-shaped envelope `kᴿ(θ) ∝ √|cos θ|` in `(kx, kz)` space, inside of which Rossby-like waves dominate over turbulence; energy accumulates at the pinch axis (the `kz` axis here, the `ky` axis in classical β-plane).

If something else in the scripts is unclear, the docstring at the top of each script generally has the math. The paper PDF and the Vallis 2017 chapter (linked above) are the most accessible references for the underlying theory.
