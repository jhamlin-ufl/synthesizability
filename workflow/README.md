# Workflow

The pipeline rebuilds every table and figure of the paper from the deposited data. Run it from
the repository root:

```bash
uv sync
uv run snakemake -c1                # all tables and Figs. 1 and 2
uv run snakemake -c1 si_figures     # Supplement fit figures (needs pdfjam)
uv run snakemake -c1 check          # compare the tables with the manuscript (PAPER_DIR)
uv run snakemake -c1 sample_readmes # rewrite samples/<formula>/README.md
```

| Output | Rule | Script |
|---|---|---|
| Table 1 (`results/tables/sample_table.tex`, `sample_phases_table.tex`) | `table_samples` | `scripts/table_samples.py` |
| Table S1 (`buttons_table.tex`) | `table_buttons` | `scripts/table_buttons.py` |
| Table S2 (`database_table.tex`) | `table_databases` | `scripts/table_databases.py` |
| Fig. 1 (`results/figures/xrd_representative.pdf`) | `fig_xrd_representative` | `scripts/fig_xrd_representative.py` |
| Fig. 2 (`results/figures/disorder_vs_stability.pdf`) | `fig_disorder_vs_stability` | `scripts/fig_disorder_vs_stability.py` |
| Supplement fit figures (`results/supplement/`) | `si_figures` | `scripts/fig_si_fits.py` |
| Hull energies (`results/hull_energies.csv`) | `hull_energies` | `scripts/hull_energies.py` |
| Fig. 2 crosses (`results/phonon_flags.csv`) | `phonon_flags` | `scripts/phonon_flags.py` |

## Folders

| Path | Content |
|---|---|
| `Snakefile` | the pipeline |
| `scripts/` | the scripts above; `refine.py`, the GSAS-II refinement driver named in every `command.txt`; `check_fit_commands.py`, which checks the paths in every `command.txt`; `sample_readmes.py`, which writes the sample READMEs; `searchmatch.py` and `searchmatch_batch.py`, the automated search-match; `chain_TiPt5Si_ZrPt5Si/`, the scripts behind two main fits (see its README); `build_snapshot.py`, which made `databases/snapshot/` from the lab repository |
| `src/synthesizability/` | library code: hull energies (`hull.py`), file parsers, plotting helpers, and the search-match (`phaseid/`, see its README) |
| `figure_inputs/` | inputs of Fig. 1: `xrd_representative_fits.csv` (one row per panel, naming its fit), `structure_insets/` (pre-rendered drawings of the refined structures), `predicted/` (the OQMD structure simulated in panel f) |

## Optional extras

```bash
uv sync --extra rietveld   # Python packages GSAS-II needs (install GSAS-II itself as below)
uv sync --extra disorder   # P_disorder classifier dependencies (CPU torch)
uv sync --extra phonons    # MatterSim + phonopy (CPU)
uv sync --extra refresh    # Materials Project client, for refreshing the database snapshot
```

`torch` comes from the PyTorch CPU index, so no CUDA libraries are installed.

## GSAS-II

Rerunning a fit needs GSAS-II, installed separately:

- conda: `conda install -c conda-forge gsas2pkg`, or
- from source: clone https://github.com/AdvancedPhotonSource/GSAS-II, build its Fortran
  extensions with meson (needs gfortran), and put the clone on `PYTHONPATH`.

The fits were run with GSAS-II built from source at commit b2e22ecc. If GSAS-II is not importable
as a package, set `GSAS2_DIR` to its `GSAS-II` source folder.
