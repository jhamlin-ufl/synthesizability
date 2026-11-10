# Sc3Pd2Cu

Table 1 row 13, candidate list A. Discussed in the Supplement, Sec. S7, under Sc3Pd2Cu.

Outcome: B2 Sc(Cu,Pd) after anneal, Sc ordered; Cu/Pd layer order of the on-hull structure absent

Category: Disordered parent of on-hull structure (firm). Not previously reported.

P_disorder = 0.88. Database structures below the hull that are dynamically stable (MatterSim): 0 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0473 | target | 0.21 | `0473_HM_Sc3CuPd2.txt` (Siemens D500, as cast); `0473_HM_Sc3CuPd2.xy` (PANalytical X'Pert MPD, as cast); `0473_HM_Sc3CuPd2_annealed.txt` (Siemens D500, annealed 1000 C 1 week) |

Scan used: `buttons/0473/0473_HM_Sc3CuPd2_annealed.txt` (Siemens D500, powder, PMMA slide (inferred), annealed 1000 C for one week).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0473_FoR_B2` | main fit | B2 Sc(Cu,Pd), Sc order eta 0.75, annealed | 4.26 |
| `0473_pred_P4mmm` | contrast | predicted P4/mmm | 5.71 |

## Database entries

Chemical system Cu-Pd-Sc. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Cu-Pd-Sc.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/Sc3CuPd2_1374854.cif`
- search-match structures: `../../databases/searchmatch/<source>/Cu-Pd-Sc/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
