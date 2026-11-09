# ScCrC2

Table 1 row 6, candidate list C. Discussed in the Supplement, Sec. S7, under ScCrC2.

Outcome: On-hull P6_3/mmc ScCrC2, nearly single phase (about 95 wt%); previously known

Category: On-hull structure (firm). Structure previously reported at this composition.

P_disorder = 0.05. Database structures below the hull that are dynamically stable (MatterSim): 2 of 2.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0512 | target | 3.00 | `0512_HM_C2CrSc.txt` (Siemens D500, as cast) |
| 0605 | remake | not recorded | `0605_HM_C2CrSc.txt` (Siemens D500, as cast) |

Scan of record: `buttons/0605/0605_HM_C2CrSc.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0605_FoR` | fit of record | P6_3/mmc ScCrC2 (Poettgen 1995) + graphite | 6.21 |
| `0605_ScCr_mixed` | contrast | Sc/Cr-mixed cell (same cell) | 7.50 |

## Database entries

Chemical system C-Cr-Sc. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/C-Cr-Sc.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/C2CrSc_14516.cif`
- search-match structures: `../../databases/searchmatch/<source>/C-Cr-Sc/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
