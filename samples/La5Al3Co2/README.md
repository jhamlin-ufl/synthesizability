# La5Al3Co2

Table 1 row 21, candidate list B. Discussed in the Supplement, Sec. S7, under La5Al3Co2.

Outcome: Known La4Co3Al3 majority; on-hull structure not identified, minority up to ~10 wt% not excluded

Category: Multiphase, on-hull structure not formed (probable). Not previously reported.

P_disorder = 0.75. Database structures below the hull that are dynamically stable (MatterSim): 0 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0483 | target | 0.03 | `0483_HM_La5Al3Co2.txt` (Siemens D500, as cast); `0483_HM_La5Al3Cu2.xy` (PANalytical X'Pert MPD, as cast) |

Scan of record: `buttons/0483/0483_HM_La5Al3Co2.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0483_FoR` | fit of record | best match: La4Co3Al3 (Tursina 2023) | 6.60 |
| `0483_pred_alone` | contrast | predicted La5Al3Co2 alone | 14.43 |

## Database entries

Chemical system Al-Co-La. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Al-Co-La.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/La5Al3Co2_1782939.cif`
- search-match structures: `../../databases/searchmatch/<source>/Al-Co-La/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
