# YSnCu

Table 1 row 5, candidate list C. Discussed in the Supplement, Sec. S7, under YSnCu.

Outcome: On-hull LiGaGe-type YSnCu formed; previously known compound

Category: On-hull structure (firm). Structure previously reported at this composition.

P_disorder = 0.15. Database structures below the hull that are dynamically stable (MatterSim): 7 of 8.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0507 | target | 0.43 | `0507_HM_CuSnY_circular.txt` (Siemens D500, as cast) |

Scan of record: `buttons/0507/0507_HM_CuSnY_circular.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0507_FoR` | fit of record | LiGaGe-type YCuSn (ordered prediction) + Y3Cu4Sn4 | 4.86 |
| `0507_CuSn_mixed_c` | contrast | Cu/Sn-mixed parent (same cell) | 8.69 |

## Database entries

Chemical system Cu-Sn-Y. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Cu-Sn-Y.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/CuSnY_1468383.cif`
- search-match structures: `../../databases/searchmatch/<source>/Cu-Sn-Y/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
