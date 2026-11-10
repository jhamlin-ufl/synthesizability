# GdFeGe3

Table 1 row 29, candidate list B. Discussed in the Supplement, Sec. S7, under GdFeGe3.

Outcome: Multiphase, three known phases; on-hull structure absent

Category: Multiphase, on-hull structure not formed (firm). Not previously reported.

P_disorder = 0.01. Database structures below the hull that are dynamically stable (MatterSim): 0 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0495 | earlier attempt | 2.26 | `0495_HM_GdFeGe3_circular.txt` (Siemens D500, as cast) |
| 0553 | target | 2.34 | `0553_HM_GdFeGe3.txt` (Siemens D500, as cast) |

Scan used: `buttons/0553/0553_HM_GdFeGe3.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0553_FoR` | main fit | GdFe0.37Ge2 + Gd0.5Fe3Ge3 + Ge | 7.52 |
| `0553_pred_alone` | contrast | predicted GdFeGe3 alone | 39.94 |

## Database entries

Chemical system Fe-Gd-Ge. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Fe-Gd-Ge.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/GdFeGe3_1364532.cif`
- search-match structures: `../../databases/searchmatch/<source>/Fe-Gd-Ge/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
