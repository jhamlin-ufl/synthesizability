# YAl3Cu2

Table 1 row 2, candidate list C. Discussed in the Supplement, Sec. S7, under YAl3Cu2.

Outcome: Predicted Al/Cu order on the CaCu5 lattice (order parameter 0.7 after anneal); other phases not identified

Category: On-hull structure (probable). Ordered member of a reported series whose known member is partly disordered.

P_disorder = 0.18. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0506 | earlier attempt | 14.31 | `0506_HM_Al3Cu2Y_circular.txt` (Siemens D500, as cast) |
| 0540 | target | 2.99 | `0540_HM_Al3Cu2Y_circular.txt` (Siemens D500, as cast) |
| 0615 | remake | 2.09 | `0615_HM_Al3Cu2Y_annealed_2x_ZB.txt` (Siemens D500, annealed 1000 C 1 week, then 500 C 12 h); `0615_HM_Al3Cu2Y_annealed_ZB.txt` (Siemens D500, annealed 1000 C 1 week) |

Scan used: `buttons/0615/0615_HM_Al3Cu2Y_annealed_ZB.txt` (Siemens D500, powder, zero-background Si plate, annealed 1000 C for one week).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0615_FoR_order` | main fit | CaCu5-type YAl3Cu2 (Al/Cu order eta 0.71) + La3Al11-type Y3(Al,Cu)11 | 11.91 |
| `0615_random` | contrast | Al/Cu random CaCu5 cell | 14.49 |

## Database entries

Chemical system Al-Cu-Y. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Al-Cu-Y.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/Al3Cu2Y_1370163.cif`
- search-match structures: `../../databases/searchmatch/<source>/Al-Cu-Y/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
