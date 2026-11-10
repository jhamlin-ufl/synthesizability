# ScPd2Cu

Table 1 row 20, candidate list A. Discussed in the Supplement, Sec. S7, under ScPd2Cu.

Outcome: fcc-derived Sc(Cu,Pd)3-type plus hexagonal phase; compositions undetermined

Category: Multiphase, on-hull structure not formed (tentative); category 3 not excluded. Not previously reported.

P_disorder = 0.75. Database structures below the hull that are dynamically stable (MatterSim): 0 of 0.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0469 | target | 1.16 | `0469_HM_ScCuPd2.txt` (Siemens D500, as cast); `0469_HM_ScCuPd2.xy` (PANalytical X'Pert MPD, as cast) |

Scan used: `buttons/0469/0469_HM_ScCuPd2.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0469_L12parent_mw5_hs` | main fit | cubic L1_2 Sc(Cu,Pd)3 + 4H hexagonal | 6.07 |
| `0469_FoR_mw5_hs` | contrast | ordered P4/mmm Sc(Cu,Pd)3 type (Alexandria, above hull) | 6.76 |
| `0469_pred_D022_w2_hs` | contrast | OQMD D0_22 prediction, DFT cell (2 % window) | 12.75 |

## Database entries

Chemical system Cu-Pd-Sc. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Cu-Pd-Sc.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/ScCuPd2_1370365.cif`
- search-match structures: `../../databases/searchmatch/<source>/Cu-Pd-Sc/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
