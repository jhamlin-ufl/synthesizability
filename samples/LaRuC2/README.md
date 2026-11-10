# LaRuC2

Table 1 row 30, candidate list C. Discussed in the Supplement, Sec. S7, under LaRuC2.

Outcome: Ru and probably La(OH)3; on-hull structure absent

Category: Multiphase, on-hull structure not formed (tentative). Not previously reported.

P_disorder = 0.00. Database structures below the hull that are dynamically stable (MatterSim): 0 of 0.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0513 | target | 6.2555 | `0513_HM_C2LaRu_circular.txt` (Siemens D500, as cast) |
| 0614 | remake | 0.62 | `0614_HM_C2LaRu.txt` (Siemens D500, as cast) |
| 0617 | remake | 2.18 | `0617_HM_C2LaRu.txt` (Siemens D500, as cast) |

Scan used: `buttons/0614/0614_HM_C2LaRu.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0614_FoR` | main fit | Ru + La(OH)3 | 7.50 |
| `0614_pred_alone` | contrast | predicted LaRuC2 alone | 40.57 |

## Database entries

Chemical system C-La-Ru. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/C-La-Ru.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/C2LaRu_1751220.cif`
- search-match structures: `../../databases/searchmatch/<source>/C-La-Ru/` (oqmd_ternary_phases, alexandria_pbe_ternary_phases, cod)
