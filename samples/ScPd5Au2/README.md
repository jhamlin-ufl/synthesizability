# ScPd5Au2

Table 1 row 14, candidate list A. Discussed in the Supplement, Sec. S7, under ScPd5Au2.

Outcome: fcc (Sc,Pd,Au) solution; Pbam order absent in range; bulk scan

Category: Disordered parent of on-hull structure (tentative). Not previously reported.

P_disorder = 0.81. Database structures below the hull that are dynamically stable (MatterSim): 0 of 0.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0477 | target | 15.57 | `0477_HM_ScPd5Au2.xy` (PANalytical X'Pert MPD, as cast) |

Scan used: `buttons/0477/0477_HM_ScPd5Au2.xy` (PANalytical X'Pert MPD, bulk piece (uncut button face), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0477_FoR_A1_c` | main fit | fcc (Sc,Pd,Au) | 15.27 |
| `0477_pred_Pbam_c` | contrast | predicted Pbam | 17.77 |

## Database entries

Chemical system Au-Pd-Sc. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Au-Pd-Sc.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/ScPd5Au2_1990792.cif`
- search-match structures: `../../databases/searchmatch/<source>/Au-Pd-Sc/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
