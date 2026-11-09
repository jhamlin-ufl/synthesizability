# Al8FeCo7

Table 1 row 9, candidate list B. Discussed in the Supplement, Sec. S7, under Al8FeCo7.

Outcome: B2 Al(Fe,Co), Al ordered (eta 0.9); Fe/Co order undetermined

Category: Disordered parent of on-hull structure (firm). Not previously reported.

P_disorder = 0.92. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0487 | earlier attempt | 9.41 | `0487_HM_Al8FeCo7.txt` (Siemens D500, as cast) |
| 0535 | target | 0.56 | `0535_HM_Al8FeCo7_circular.txt` (Siemens D500, as cast) |

Scan of record: `buttons/0535/0535_HM_Al8FeCo7_circular.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0535_FoR_B2` | fit of record | B2 Al(Fe,Co), Fe/Co mixed, Al order eta 0.9 | 5.91 |
| `0535_pred_Al_eta0.90` | contrast | predicted Fe/Co-ordered Pm-3m, same Al order | 5.91 |

## Database entries

Chemical system Al-Co-Fe. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Al-Co-Fe.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/Al8FeCo7_2154926.cif`
- search-match structures: `../../databases/searchmatch/<source>/Al-Co-Fe/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
