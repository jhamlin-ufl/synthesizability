# Automated search-match

Code of the automated phase identification described in Supplement Section S2. For one button it
builds a library of candidate structures for the ternary system, ranks combinations of one to
three phases against the measured pattern, and refines the best-ranked combinations in GSAS-II.

| Module | Role |
|---|---|
| `library.py` | builds the library: OQMD, Materials Project and Alexandria entries of the system and its subsystems, COD structures, the predicted structures, a fully disordered version of each ordered structure, partially disordered versions (MatterSim-relaxed), Vegard solid solutions between isostructural known binaries, and random fcc/bcc/hcp cells at the nominal composition |
| `mixtures.py` | enumerates single phases, pairs and triples whose compositions can combine to the nominal one |
| `searchmatch.py` | scales each calculated pattern by up to ±2 % in cell, fits combinations by non-negative least squares and scores them |
| `refine.py` | GSAS-II refinement of the best-ranked combinations |
| `structure_types.py` | groups entries of the same structure type |

## Inputs

All in `databases/searchmatch/` (set `SYNTH_SEARCH_DIR` to use another location), for the chemical
systems of the 34 candidates only; see `databases/README.md`. Measured patterns come from
`samples/<formula>/buttons/<NNNN>/`, and the nominal formula of a button from `tables/buttons.csv`.

ICDD PDF cards are not included. In the study, the library also held the card CIFs exported
from our JADE analysis. Users with a PDF licence can add their own: set `SYNTH_PDF_CARD_DIR` to a
folder of card CIFs, and every card whose elements lie in the system is added. Without it, the
library is built from the open sources alone.

## Running

```bash
uv sync --extra rietveld                      # GSAS-II itself: see workflow/README.md
uv run python workflow/scripts/searchmatch.py --sample 561 --scan 0561_HM_ZrPt5Si_annealed_ZB.txt --no-refine
uv run python workflow/scripts/searchmatch.py --sample 561 --scan 0561_HM_ZrPt5Si_annealed_ZB.txt --top 15 --diverse
uv run python workflow/scripts/searchmatch_batch.py --name test --samples 561 493 --top 15 --diverse
```

Results go to `results/searchmatch/<button>/`; the library cache to `results/searchmatch/library/`.
`--no-refine` (ranking only) needs no GSAS-II. The study used `--top 15 --diverse`.

The search proposes candidates; the final fits in `samples/*/fits/` were then chosen and refined
as described in Supplement Sections S2 and S3, each with its own `command.txt`.
