# La2PdAu2

Table 1 row 32, candidate list A. Discussed in the Supplement, Sec. S7, under La2PdAu2.

Outcome: Not the on-hull structure; multiphase, no phase identified

Category: On-hull structure did not form but structure(s) unresolved (firm). Not previously reported.

P_disorder = 0.88. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0475 | target | 0.9 | `0475_HM_La2PdAu2.txt` (Siemens D500, as cast); `0475_HM_La2PdAu2.xy` (PANalytical X'Pert MPD, as cast) |

Scan of record: `buttons/0475/0475_HM_La2PdAu2.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0475_best_hs` | best fit found (not a satisfactory description; not shown in the Supplement) | best three-phase description (not a phase assignment) | 14.31 |
| `0475_pred_alone` | contrast | predicted La2PdAu2 alone | 16.90 |

## Database entries

Chemical system Au-La-Pd. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Au-La-Pd.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/La2PdAu2_1980865.cif`
- search-match structures: `../../databases/searchmatch/<source>/Au-La-Pd/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
