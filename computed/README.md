# Computed inputs

Model results the pipeline reads as they are. They were computed in the group's lab repository,
because they need extra software or model weights that are not redistributed here.

| Path | What |
|---|---|
| `p_disorder.csv` | P_disorder of the 34 targets |
| `phonons/` | dynamical stability of the database structures near the hull |

The hull energies are not stored here. The pipeline computes them from `databases/snapshot/`
(`workflow/scripts/hull_energies.py`, written to `results/hull_energies.csv`).

## P_disorder

Probability of substitutional disorder from the classifier of Jakob et al. (2025), a recurrent
network that uses the composition alone (their released `hyperopt_best_model.pt` and
`hyperopt_config.json`). To recompute it, install their code with `elementembeddings==0.6.1` and
torch (CPU is enough; `uv sync --extra disorder`) and score the 34 formulas. The 0.6.2 wheel of
`elementembeddings` lacks its data folder, and every formula then fails without an error. Check
values: Al3Cu2Y 0.1809, CuSnY 0.1502, Hf2MoIr 0.8958.

## phonons/

PhononBench protocol (Han et al. 2026, arXiv:2512.21227) with the MatterSim-v1.0.0-5M potential
(mattersim 1.2.5, phonopy 4.8.0, seekpath 2.2.2, ase 3.29.0, torch 2.14.1 CPU). Each structure is
relaxed with its space group fixed (FIRE, fmax 0.005 eV/A). Force constants come from a 2x2x2
supercell with 0.01 A displacements. A structure is unstable if any frequency on the seekpath band
path is below -0.1 THz.

| File | Content |
|---|---|
| `candidates.csv` | every OQMD, MP and Alexandria PBE entry at a target composition within 50 meV/atom of the hull, grouped into distinct structures (134 entries, 84 structures) |
| `structures/`, `relaxed_5M/` | the input cells and the MatterSim-relaxed cells |
| `phonon_results_5M.csv` | one row per structure: lowest band and mesh frequency, verdict (the paper's set) |
| `phonon_results_1M.csv` | the same with the 1M checkpoint (sensitivity check; 10 entries disagree) |
| `summary_5M_paper.csv`, `polymorphs_5M.csv` | summaries per target and per structure |
| `dynamic_flags_5M_paper.csv` | the Fig. 2 crosses; `results/phonon_flags.csv` must equal it |

The run takes about 40 CPU-minutes (`uv sync --extra phonons`). Verdicts with a band minimum near
-0.1 THz can change with the MatterSim version.

## Code

The scripts that produced these files are in the lab repository and are not part of this
release: `compute_disorder_probabilities.py` and `src/synthesizability/disorder.py` for
P_disorder; `enumerate_phonon_candidates.py`, `run_phonon_stability.py` and
`phonon_stability_core.py` for the phonons. In this release `workflow/scripts/phonon_flags.py`
recomputes the hull energies of the phonon candidates from `databases/snapshot/` and checks the
Fig. 2 crosses.
