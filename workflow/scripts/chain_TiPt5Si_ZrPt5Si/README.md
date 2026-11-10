# Refinement chain of the TiPt5Si (0493) and ZrPt5Si (0561) main fits

These two main fits were made before the general driver `workflow/scripts/refine.py` existed, by a
three-step chain run on the workstation:

1. `refine_0493_annealed_tisipt5.py` — the predicted P4/mmm phase alone, from its OQMD cell;
2. `refine_0493_annealed_two_phase.py` — adds the second phase (fcc for 0493, D0_24 ZrPt3 for 0561);
3. `refine_0493_annealed_polish.py --sample 493|561` — one block at a time, each kept only if it lowers
   R_wp by 0.05 and stays physical: axial-divergence asymmetry, sample transparency, Pt Uiso,
   surface roughness.

The scripts write each step to `results/refits/chain/<sample>_<step>/` and read the starting structures
from the search-match library (`results/searchmatch/library/<space>/`, built by `workflow/scripts/searchmatch.py`). The resulting
GSAS-II projects (`samples/TiPt5Si/fits/0493_annealed_polish/work/fit.gpx`,
`samples/ZrPt5Si/fits/0561_annealed_polish/work/fit.gpx`) are deposited and re-refine in GSAS-II.

For both samples the Supplement also shows an ordered fit made with `workflow/scripts/refine.py`
(`.../ordered_matched/`, with its `command.txt`), which comes within 0.4 percentage points of the
main fit in R_wp and is the one compared with the site-mixed model.
