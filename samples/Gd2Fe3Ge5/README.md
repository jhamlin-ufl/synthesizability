# Gd2Fe3Ge5

Table 1 row 28, candidate list B. Discussed in the Supplement, Sec. S7, under Gd2Fe3Ge5.

Outcome: Multiphase, three known phases; on-hull structure absent

Category: Multiphase, on-hull structure not formed (firm). Not previously reported.

P_disorder = 0.01. Database structures below the hull that are dynamically stable (MatterSim): 2 of 2.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0496 | earlier attempt | 1.83 | `0496_HM_Gd2Fe3Ge5_circular.txt` (Siemens D500, as cast) |
| 0554 | target | 4.49 | `0554_HM_Gd2Fe3Ge5.txt` (Siemens D500, as cast) |
| 0557 | remake | 4.96 | none |

Scan used: `buttons/0554/0554_HM_Gd2Fe3Ge5.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0554_tex_FoR` | main fit | GdFe0.37Ge2 + GdFe2Ge2 (preferred orientation) + Gd0.5Fe3Ge3 + Ge | 9.40 |
| `0554_pred_alone` | contrast | predicted Gd2Fe3Ge5 alone | 38.49 |

## Database entries

Chemical system Fe-Gd-Ge. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Fe-Gd-Ge.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/Gd2Fe3Ge5_1800260.cif`
- search-match structures: `../../databases/searchmatch/<source>/Fe-Gd-Ge/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
