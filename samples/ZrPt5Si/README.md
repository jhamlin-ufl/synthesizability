# ZrPt5Si

Table 1 row 1, candidate list B. Discussed in the Supplement, Sec. S7, under ZrPt5Si.

Outcome: On-hull ZrPt5Si, 97 wt% after anneal (ZrPt3 3 wt%); new compound

Category: On-hull structure (probable). Not previously reported.

P_disorder = 0.20. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0482 | earlier attempt | 0.27 | `0482_HM_ZrSiPt5.xy` (PANalytical X'Pert MPD, as cast) |
| 0561 | target | 0.99 | `0561_HM_ZrPt5Si_annealed_ZB.txt` (Siemens D500, annealed 1000 C 1 week); `0561_HM_ZrSiPt5.txt` (Siemens D500, as cast) |

Scan of record: `buttons/0561/0561_HM_ZrPt5Si_annealed_ZB.txt` (Siemens D500, powder, zero-background Si plate, annealed 1000 C for one week).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0561_annealed_polish` | fit of record | ZrPt5Si annealed, ordered P4/mmm | 8.77 |
| `0561_ordered_matched` | contrast | ZrPt5Si annealed, ordered P4/mmm, refined the same way as the mixed model | 8.43 |
| `0561_sitemixed` | contrast | ZrPt5Si annealed, Zr/Si/Pt site-mixed parent | 23.65 |
| `0561_ordered_noZrPt3` | contrast | ZrPt5Si annealed, ordered P4/mmm without ZrPt3, refined the same way | 10.45 |
| `0561_ascast_joint_zero` | as-cast scan | ZrPt5Si as cast | 7.72 |
| `0561_exchange_md` | quoted in the Supplement text, not plotted | ordered ZrPt5Si with refined Zr/Si exchange (about 1 %), with texture | 8.42 |
| `0561_exchange_nomd` | quoted in the Supplement text, not plotted | ordered ZrPt5Si with refined Zr/Si exchange, without texture | 11.08 |

## Database entries

Chemical system Pt-Si-Zr. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Pt-Si-Zr.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/ZrSiPt5_2146113.cif`
- search-match structures: `../../databases/searchmatch/<source>/Pt-Si-Zr/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
