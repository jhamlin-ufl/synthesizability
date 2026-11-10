# Samples

One folder per target, named by its formula as printed in Table 1. Each folder has a
`README.md` with the outcome, the buttons, the scans and the fits of that target.

```
samples/<formula>/
  README.md
  buttons/<NNNN>/          one arc-melted button
  fits/<NNNN>_<fit>/       one GSAS-II refinement of a scan of button NNNN
```

A target can have several buttons: earlier attempts, remakes and second buttons
(`tables/buttons.csv`). The scan used is not always from the first button. For example, the
YAl3Cu2 scan used is from the remake 0615.

## Buttons

Each `buttons/<NNNN>/` folder holds the synthesis record and the measured patterns.

| File | Content |
|---|---|
| `SYNTHESIS` or `SYNTHESIS.txt` | the synthesis record: starting materials, target and weighed masses, mass loss |
| `STATUS` | short notes, where kept: superconductivity check, XRD, candidate list (older list names, see `tables/lists/README.md`) |
| `*.txt` | Siemens D500 scan, Bruker RAW4 text export with its header (Cu K-alpha) |
| `*.xy` | PANalytical X'Pert MPD scan, two columns: 2-theta (deg) and counts |

The scan file names keep the formula spelling of the lab records (for example `Cu5GePd6` for
Pd6Cu5Ge). `tables/scans.csv` gives the instrument, mount, range and state of every scan.

## Fits

Each `fits/<NNNN>_<fit>/` folder is one refinement. `tables/si_fits.csv` says which fits the
Supplement shows and what each one tests. The sample READMEs list the others.

| File | Content |
|---|---|
| `fit.txt` | columns `two_theta obs calc background diff` on the refined range |
| `ticks.json` | Bragg positions (deg) of each phase, in refinement order |
| `summary.json` | `record.wR` (R_wp, %), `record.state.wt` (weight fractions, in `ticks.json` order), cells and atoms; `stages` lists each refinement step |
| `scan.xy` | the scan as given to GSAS-II |
| `fit.gpx` or `work/fit.gpx` | the GSAS-II project |
| `start/` | the starting CIFs |
| `command.txt` | the command that made the fit |

`<fit>` names: `FoR` main fit, `pred` the predicted structure, `_hs` sample displacement held
at zero because the refined value was unphysical. Phase names give the source of each starting
model: `oqmd`, `mp`, `alex` (Alexandria), `cod`, a literature reference (for example
`Poettgen1995`), `pred` (the predicted cell), `_dis` (its site-disordered parent), `pm` (a solid
solution between two end members). Every starting structure comes from an open source or was
entered from a publication (provenance in the CIF header). None comes from an ICDD card.

### Rerunning a fit

Install GSAS-II (see `workflow/README.md`), run `uv sync --extra rietveld`, then run the command
in `command.txt` from the repository root. The new fit goes to `results/refits/`.
`python workflow/scripts/check_fit_commands.py` checks that every file named in every
`command.txt` exists.

Three fits were not made with `workflow/scripts/refine.py`. `TiPt5Si/fits/0493_annealed_polish`
and `ZrPt5Si/fits/0561_annealed_polish` came from the scripts in
`workflow/scripts/chain_TiPt5Si_ZrPt5Si/`. `GdPt2Ge/fits/0480_ascast_90321` came from the
automated search-match, and its `summary.json` uses the older keys `rec.wR` and `rec.phases[].wt`.
Their GSAS-II projects re-refine in GSAS-II.

GSAS-II stores absolute file paths inside each `.gpx`. They point to the computers the fits were
run on. They do not stop a project from opening, because the scan and phases are stored in it.

## Files left out

`EXCLUDED.md` lists the files of the lab records that are not deposited: ICDD card exports, MDI
JADE work files, intermediate structure files and a few others.
