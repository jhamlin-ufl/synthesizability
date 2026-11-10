# Hf2MoIr

Table 1 row 10, candidate list A. Discussed in the Supplement, Sec. S7, under Hf2MoIr.

Outcome: bcc Hf-Mo-Ir after anneal, partial Hf order; Mo/Ir order absent; Ti2Ni-type Hf2Ir 11 wt%

Category: Disordered parent of on-hull structure (probable). Not previously reported.

P_disorder = 0.90. Database structures below the hull that are dynamically stable (MatterSim): 2 of 2.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0466 | target | 0.55 | `0466_HM_Hf2MoIr 1.xy` (PANalytical X'Pert MPD, as cast); `0466_HM_Hf2MoIr_annealed.txt` (Siemens D500, annealed 1000 C 1 week); `0466_HM_Hf2MoIr_annealed_ZB.txt` (Siemens D500, annealed 1000 C 1 week) |

Scan used: `buttons/0466/0466_HM_Hf2MoIr_annealed_ZB.txt` (Siemens D500, powder, zero-background Si plate, annealed 1000 C for one week).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0466_FoR_B2` | main fit | bcc, partial Hf order (eta 0.64) + Hf2Ir | 7.31 |
| `0466_pred_L21` | contrast | predicted L2_1 Hf2MoIr + Hf2Ir | 16.89 |

## Database entries

Chemical system Hf-Ir-Mo. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Hf-Ir-Mo.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/Hf2MoIr_1013867.cif`
- search-match structures: `../../databases/searchmatch/<source>/Hf-Ir-Mo/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
