# synthesizability

Data, fits and code for

> G. D. Di Gianluca *et al.*, *Experimental Evaluation of High-Throughput Structure Predictions*,
> [arXiv:2611.09671](https://arxiv.org/abs/2611.09671) (2026).

We arc-melted 34 ternary compounds predicted to be stable or nearly stable and identified the
phases that formed by powder X-ray diffraction and Rietveld refinement.

## Contents

| Path | What it holds |
|---|---|
| `Supplement.pdf` | the Supplementary Information |
| `samples/` | one folder per target, named by its formula in Table 1: synthesis records, X-ray scans and GSAS-II fits |
| `tables/` | the targets, buttons, scans and fits as CSV tables; the candidate lists |
| `computed/` | P_disorder and the dynamical-stability (phonon) results |
| `databases/` | OQMD, Materials Project and Alexandria entries for each chemical system, and the search-match structures |
| `instruments/` | GSAS-II instrument profiles of the two diffractometers |
| `results/` | the generated tables and figures |
| `workflow/` | the pipeline (`Snakefile`), its scripts and library code |

## Finding a sample

Open `samples/<formula>/README.md`, for example `samples/ZrPt5Si/README.md`. It gives the Table 1
row and outcome, the buttons and scans, the main fit and the contrast fits, and where the
database entries for the chemical system are. `samples/README.md` explains the file formats and
how to rerun a fit.

## Rebuilding the tables and figures

```bash
uv sync
uv run snakemake -c1              # Table 1, Tables S1 and S2, Figs. 1 and 2
uv run snakemake -c1 si_figures   # the Supplement's fit figures (needs pdfjam)
```

The outputs go to `results/`. `workflow/README.md` lists which script makes which output.

## Not included

ICDD PDF-5+ cards and MDI JADE work files are licensed and not redistributed. No deposited fit
uses an ICDD card. `samples/EXCLUDED.md` lists the files left out.

## Licence and citation

Code (`workflow/`): MIT (`LICENSE`). Data, figures and text: CC BY 4.0 (`LICENSE-DATA`). The
database extracts keep the licences of their sources (see `databases/README.md`).

Please cite the paper and this archive, DOI
[10.5281/zenodo.19685017](https://doi.org/10.5281/zenodo.19685017) (see `CITATION.cff`).
