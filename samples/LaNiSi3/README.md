# LaNiSi3

Table 1 row 24, candidate list C. Discussed in the Supplement, Sec. S7, under LaNiSi3.

Outcome: Multiphase, known compounds; on-hull structure absent

Category: Multiphase, on-hull structure not formed (firm). Not previously reported.

P_disorder = 0.07. Database structures below the hull that are dynamically stable (MatterSim): 0 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0514 | target | 2.2 | `0514_HM_LaNiSi3_circular.txt` (Siemens D500, as cast) |

Scan used: `buttons/0514/0514_HM_LaNiSi3_circular.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0514_FoR` | main fit | La3Ni3Si7 + LaNiSi2 + Si + LaNi6Si6 (Pani 2014 isotype model) | 10.64 |
| `0514_pred_alone` | contrast | predicted LaNiSi3 alone | 41.27 |

## Database entries

Chemical system La-Ni-Si. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/La-Ni-Si.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/LaNiSi3_1364676.cif`
- search-match structures: `../../databases/searchmatch/<source>/La-Ni-Si/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
