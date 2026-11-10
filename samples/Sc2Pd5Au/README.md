# Sc2Pd5Au

Table 1 row 16, candidate list A. Discussed in the Supplement, Sec. S7, under Sc2Pd5Au.

Outcome: L1_2 Sc(Pd,Au)3; Pbam Pd/Au order absent in range

Category: Disordered parent of on-hull structure (probable). Not previously reported.

P_disorder = 0.60. Database structures below the hull that are dynamically stable (MatterSim): 3 of 3.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0478 | target | 0.46 | `0478_HM_Sc2Pd5Au.xy` (PANalytical X'Pert MPD, as cast) |
| 0610 | remake | 1.01 | `0610_HM_Sc2Pd5Au.txt` (Siemens D500, as cast) |

Scan used: `buttons/0610/0610_HM_Sc2Pd5Au.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0610_FoR_L12_c` | main fit | L1_2 Sc(Pd,Au)3, Sc order eta 0.9 | 7.16 |
| `0610_pred_Pbam` | contrast | predicted Pbam | 8.98 |

## Database entries

Chemical system Au-Pd-Sc. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Au-Pd-Sc.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/Sc2Pd5Au_1992186.cif`
- search-match structures: `../../databases/searchmatch/<source>/Au-Pd-Sc/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
