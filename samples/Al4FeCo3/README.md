# Al4FeCo3

Table 1 row 11, candidate list B. Discussed in the Supplement, Sec. S7, under Al4FeCo3.

Outcome: B2 Al(Fe,Co), Al ordered (eta 0.9); Fe/Co order undetermined

Category: Disordered parent of on-hull structure (firm). Not previously reported.

P_disorder = 0.89. Database structures below the hull that are dynamically stable (MatterSim): 3 of 3.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0486 | earlier attempt | 3.94 | `0486_HM_Al4FeCo3.txt` (Siemens D500, as cast) |
| 0534 | remake | 3.33 | `0534_HM_Al4FeCo3_circular.txt` (Siemens D500, as cast) |
| 0544 | remake | 10.88 | none |
| 0545 | target | 0.82 | `0545_HM_Al4FeCo3_circular.txt` (Siemens D500, as cast) |

Scan of record: `buttons/0545/0545_HM_Al4FeCo3_circular.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0545_FoR_B2_c` | fit of record | B2 Al(Fe,Co), Fe/Co mixed, Al order eta 0.9 | 5.78 |
| `0545_pred_Al_eta0.90` | contrast | predicted Fe/Co-ordered Im-3m, same Al order | 5.78 |

## Database entries

Chemical system Al-Co-Fe. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Al-Co-Fe.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/Al4FeCo3_2146849.cif`
- search-match structures: `../../databases/searchmatch/<source>/Al-Co-Fe/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
