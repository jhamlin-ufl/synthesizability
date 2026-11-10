# GdFeGe2

Table 1 row 27, candidate list B. Discussed in the Supplement, Sec. S7, under GdFeGe2.

Outcome: Two known phases; on-hull structure absent

Category: Multiphase, on-hull structure not formed (firm). Not previously reported.

P_disorder = 0.01. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0497 | earlier attempt | 5.40 | `0497_HM_GdFeGe2_circular.txt` (Siemens D500, as cast) |
| 0555 | target | 1.49 | `0555_HM_GdFeGe2.txt` (Siemens D500, as cast) |

Scan used: `buttons/0555/0555_HM_GdFeGe2.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0555_tex_FoR` | main fit | GdFe0.37Ge2 + GdFe2Ge2 (preferred orientation) + Gd0.5Fe3Ge3 + Gd11Ge10 | 7.26 |
| `0555_pred_alone` | contrast | predicted GdFeGe2 alone | 37.39 |

## Database entries

Chemical system Fe-Gd-Ge. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Fe-Gd-Ge.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/GdFeGe2_1377698.cif`
- search-match structures: `../../databases/searchmatch/<source>/Fe-Gd-Ge/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
