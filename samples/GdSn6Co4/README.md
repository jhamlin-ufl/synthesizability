# GdSn6Co4

Table 1 row 25, candidate list B. Discussed in the Supplement, Sec. S7, under GdSn6Co4.

Outcome: Multiphase, known compounds, GdCo6Sn6 majority; on-hull structure absent

Category: Multiphase, on-hull structure not formed (firm). Not previously reported.

P_disorder = 0.04. Database structures below the hull that are dynamically stable (MatterSim): 0 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0490 | earlier attempt | 11.09 | `0490_HM_GdCo4Sn6_circular.txt` (Siemens D500, as cast) |
| 0537 | target | 0.5 | `0537_HM_GdCo4Sn6_circular.txt` (Siemens D500, as cast) |

Scan used: `buttons/0537/0537_HM_GdCo4Sn6_circular.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0537_FoR` | main fit | GdCo6Sn6-type + Gd3Co4Sn13 + GdCo_xSn_2 + CoSn + Sn | 7.77 |
| `0537_pred_alone` | contrast | predicted GdSn6Co4 alone | 36.47 |

## Database entries

Chemical system Co-Gd-Sn. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Co-Gd-Sn.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/GdCo4Sn6_1368900.cif`
- search-match structures: `../../databases/searchmatch/<source>/Co-Gd-Sn/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
