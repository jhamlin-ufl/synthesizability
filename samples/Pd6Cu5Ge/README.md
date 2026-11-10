# Pd6Cu5Ge

Table 1 row 31, candidate list A. Discussed in the Supplement, Sec. S7, under Pd6Cu5Ge.

Outcome: Not Cmmm; B2-type cut face, broad fcc after grinding; unresolved

Category: On-hull structure did not form but structure(s) unresolved (tentative). Not previously reported.

P_disorder = 0.93. Database structures below the hull that are dynamically stable (MatterSim): 0 of 2.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0465 | target | 0.30 | `0465_HM_Cu5GePd6.xy` (PANalytical X'Pert MPD, as cast) |
| 0568 | remake | 1.61 | `0568_HM_Cu5GePd6.txt` (Siemens D500, as cast) |
| 0608 | remake | 1.08 | `0608_HM_Cu5GePd6.txt` (Siemens D500, as cast) |

Scan used: `buttons/0568/0568_HM_Cu5GePd6.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0568_shift_best_hs` | best fit found (not a satisfactory description; not shown in the Supplement) | best fit: broad fcc + bcc solid solutions | 8.72 |
| `0568_shift_pred_hs` | contrast | predicted Cmmm alone | 14.39 |

## Database entries

Chemical system Cu-Ge-Pd. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Cu-Ge-Pd.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/Cu5GePd6_1992694.cif`
- search-match structures: `../../databases/searchmatch/<source>/Cu-Ge-Pd/` (oqmd_ternary_phases, alexandria_pbe_ternary_phases, cod)
