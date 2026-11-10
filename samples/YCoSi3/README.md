# YCoSi3

Table 1 row 8, candidate list C. Discussed in the Supplement, Sec. S7, under YCoSi3.

Outcome: Cmmm polymorph (+8 meV/atom), nearly single phase after anneal; on-hull I4mm structure absent

Category: Above-hull polymorph (probable). Composition previously reported, without a structure.

P_disorder = 0.01. Database structures below the hull that are dynamically stable (MatterSim): 2 of 3.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0504 | earlier attempt | 2.94 | `0504_HM_CoSi3Y_circular.txt` (Siemens D500, as cast) |
| 0539 | target | 0.46 | `0539_HM_CoSi3Y_annealed.txt` (Siemens D500, annealed 1000 C 1 week); `0539_HM_CoSi3Y_circular.txt` (Siemens D500, as cast); `0539_HM_YCoSi3_annealed_ZB.txt` (Siemens D500, annealed 1000 C 1 week) |

Scan used: `buttons/0539/0539_HM_YCoSi3_annealed_ZB.txt` (Siemens D500, powder, zero-background Si plate, annealed 1000 C for one week).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0539_cmmm_thin_genstrain` | main fit | Cmmm YCoSi3 (SmNiGe3 type), single phase | 12.87 |
| `0539_i4mm_oqmd_thin_genstrain` | contrast | predicted I4mm alone | 53.60 |

## Database entries

Chemical system Co-Si-Y. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Co-Si-Y.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/CoSi3Y_1475971.cif`
- search-match structures: `../../databases/searchmatch/<source>/Co-Si-Y/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
