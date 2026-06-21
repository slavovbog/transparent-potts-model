# transparent-potts-model

Julia code for exact mean-field analysis and Monte Carlo simulations of a $q$-state Potts model with many-body couplings on a complete graph.

Companion repository for [*Collective communication in a transparent world*](https://doi.org/...) (Phys. Rev. E; Akara-pipattana, Nechaev & Slavov).

For further information, contact Bogdan Slavov (slavovbog@gmail.com).

## Quick start

Requires [Julia](https://julialang.org/) ≥ 1.11.

```bash
git clone https://github.com/slavovbog/transparent-potts-model.git
cd transparent-potts-model
julia --project -e 'using Pkg; Pkg.instantiate()'
```

**Figs. 6–7** use bundled MC trajectories and run in about a minute:

```bash
julia --project scripts/make_figures_6-7.jl
```

PDFs are written to `figures/` (created automatically).

`julia --project scripts/check_representations.jl` is a sanity check for the Hamiltonian.

## Model notation

The code follows the notation in `manuscript/main.tex`. Main symbols:

| Symbol | Meaning |
|--------|---------|
| $q$ | number of colors / communities |
| $n$ | number of spins (`Inf` for the thermodynamic limit) |
| $T$ | temperature (energy units; not $\beta$) |
| $J_k$ | $k$-body coupling; `J[k-1]` in code is $J_k$ for $k \ge 2$ |
| $n_p$ | occupation number of color $p$ (`ns[p]`) |
| $c_p$ | concentration $n_p/n$ |
| $r$ | maximum interaction order ($k = 2,\ldots,r$) |

For $r=3$, couplings are `J = (J_2, J_3)`.

## Data

The Julia workflows in `src/workflows/` read and write JLD2 files through an in-memory cache (`get_or_run!` / `parallel_run!`).

### Included in the repository

| File | Description |
|------|-------------|
| `trajectories/traj_n_steps=5000.jld2` | MC concentration trajectories $c_p(t)$ for Figs. 6–7 |

Each record has the form `(input = …, output = (traj = …,))`.

- **`input`** — model parameters $(q, n, T, J)$, MC settings (`n_steps`, `burn_in`, `thin`), and `replicate_id` with deterministic `seed` (via [StableRNGs.jl](https://github.com/JuliaRandom/StableRNGs.jl)).
- **`output.traj`** — named tuple with `c` ($q \times n_{\text{recorded}}$ matrix of $c_p$), `t_sweeps`, and `accept_rate`.

Seeds are computed by `trajectory_seed(config, replicate_id)` from the model settings and replicate index. The same `replicate_id` always yields the same trajectory.

Records are keyed `"1"`, `"2"`, … in the order produced by `scripts/generate_trajectories.jl`.

### Generated locally

| Directory | Contents | Produced by |
|-----------|----------|-------------|
| `c_min_J_2_J_3/` | $c_{\min}(J_2, J_3)$ on a grid | `scripts/generate_concentration_minima.jl` |

Default grid: $J_2, J_3 \in [-200, 200]$ (401×401 = 160 801 points per $(q, n, T)$). One JLD2 file is written per parameter set, named e.g. `c_min_q=3.0_n=Inf_T=1.0.jld2`.

Each record stores `(input = …, output = (c_min = …,))` where `c_min` is a vector of minimizer locations (length 1 for a unique minimum, longer if degenerate within tolerance).


## Reproducing the figures

All figure scripts write PDFs to `figures/`.

Use `-p auto` to launch one worker per logical CPU. On a laptop, Figs. 1–5 data generation may take hours.

**Figs. 1–5** require a $401 \times 401$ grid of $c_{\min}(J_2,J_3)$ values. This is a heavy parallel job; to regenerate locally run:

```bash
julia --project -p auto scripts/generate_concentration_minima.jl
julia --project scripts/make_figures_1-5.jl
```

New **Trajectories** can be generated via:

```bash
julia --project -p auto scripts/generate_trajectories.jl
julia --project scripts/make_figures_6-7.jl
```

## Citation

If you use this code, please cite the paper:

```bibtex
@article{AkaraPipattanaNechaevSlavov2026,
  author  = {Akara-pipattana, Pawat and Nechaev, Sergei and Slavov, Bogdan},
  title   = {Collective communication in a transparent world: Phase transitions in a many-body Potts model and a social--quantum correspondence},
  journal = {Phys. Rev. E},
  year    = {2026},
  note    = {in preparation}
}
```