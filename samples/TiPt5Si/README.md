# TiPt5Si

Table 1 row 4, candidate list B. Discussed in the Supplement, Sec. S7, under TiPt5Si.

Outcome: On-hull TiPt5Si, 91 wt% after anneal; new compound

Category: On-hull structure (probable). Not previously reported.

P_disorder = 0.18. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0493 | target | 1.87 | `0493_HM_TiPt5Si_annealed_ZB.txt` (Siemens D500, annealed 1000 C 1 week); `0493_HM_TiSiPt5.txt` (Siemens D500, as cast) |
| 0618 | second button | not recorded | `0618_HM_TiPt5Si_ZB.txt` (Siemens D500, unclear (SYNTHESIS mentions a 600 C anneal of the filed powder, duration blank)) |

Scan of record: `buttons/0493/0493_HM_TiPt5Si_annealed_ZB.txt` (Siemens D500, powder, zero-background Si plate, annealed 1000 C for one week).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0493_annealed_polish` | fit of record | TiPt5Si annealed, ordered P4/mmm | 10.58 |
| `0493_ordered_matched` | contrast | TiPt5Si annealed, ordered P4/mmm, refined the same way as the mixed model | 10.57 |
| `0493_sitemixed` | contrast | TiPt5Si annealed, Ti/Si/Pt site-mixed parent | 24.55 |
| `0493_ascast_joint_zero` | as-cast scan | TiPt5Si as cast | 6.15 |

## Database entries

Chemical system Pt-Si-Ti. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Pt-Si-Ti.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/TiSiPt5_2146105.cif`
- search-match structures: `../../databases/searchmatch/<source>/Pt-Si-Ti/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
