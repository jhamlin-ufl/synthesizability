# Sc2PdCu

Table 1 row 12, candidate list A. Discussed in the Supplement, Sec. S7, under Sc2PdCu.

Outcome: B2 Sc(Cu,Pd), Sc ordered; Heusler Cu/Pd order absent

Category: Disordered parent of on-hull structure (probable). Not previously reported.

P_disorder = 0.89. Database structures below the hull that are dynamically stable (MatterSim): 0 of 2.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0474 | target | 0.09 | `0474_HM_Sc2CuPd.xy` (PANalytical X'Pert MPD, as cast); `0474_HM_Sc2CuPd_annealed.txt` (Siemens D500, annealed 1000 C 1 week) |

Scan used: `buttons/0474/0474_HM_Sc2CuPd.xy` (PANalytical X'Pert MPD, bulk piece (cut face), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0474_FoR_B2` | main fit | B2 Sc(Cu,Pd) | 8.06 |
| `0474_pred_L21` | contrast | predicted L2_1 Heusler | 10.21 |
| `0474_alex_P4mmm` | contrast | Alexandria P4/mmm (on-hull in Alexandria) | 8.65 |

## Database entries

Chemical system Cu-Pd-Sc. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Cu-Pd-Sc.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/Sc2CuPd_363099.cif`
- search-match structures: `../../databases/searchmatch/<source>/Cu-Pd-Sc/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
