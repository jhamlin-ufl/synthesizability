# ScPd2Au

Table 1 row 19, candidate list A. Discussed in the Supplement, Sec. S7, under ScPd2Au.

Outcome: Hexagonal Sc(Pd,Au)3 majority; fcc-stacked Sc(Pd,Au)3 minority is the on-hull structure or its L1_2 parent

Category: Other solid solution (probable). Not previously reported.

P_disorder = 0.56. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0476 | target | 0.18 | `0476_HM_ScPd2Au.txt` (Siemens D500, as cast); `0476_HM_ScPd2Au.xy` (PANalytical X'Pert MPD, as cast) |

Scan used: `buttons/0476/0476_HM_ScPd2Au.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0476_FoR_L12at3997_hs` | main fit | hexagonal D0_24 Sc(Pd,Au)3 + L1_2 Sc(Pd,Au)3 | 7.85 |
| `0476_pred_cubic_at3997_hs` | contrast | predicted Pd/Au order at the L1_2 cell | 8.51 |
| `0476_pred_for_L12_hs` | contrast | predicted P4/mmm, DFT cell (c/a 0.953) | 18.57 |

## Database entries

Chemical system Au-Pd-Sc. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Au-Pd-Sc.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/ScPd2Au_1988846.cif`
- search-match structures: `../../databases/searchmatch/<source>/Au-Pd-Sc/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
