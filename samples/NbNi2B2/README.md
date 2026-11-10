# NbNi2B2

Table 1 row 26, candidate list C. Discussed in the Supplement, Sec. S7, under NbNi2B2.

Outcome: Multiphase: Ni3B, NbB2, tau-boride Nb2Ni21B6; on-hull structure absent

Category: Multiphase, on-hull structure not formed (probable). Not previously reported.

P_disorder = 0.03. Database structures below the hull that are dynamically stable (MatterSim): 0 of 0.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0510 | target | 0.50 | `0510_HM_B2NbNi2_circular.txt` (Siemens D500, as cast) |

Scan used: `buttons/0510/0510_HM_B2NbNi2_circular.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0510_FoR_c` | main fit | Ni3B + NbB2 + tau-boride Nb2Ni21B6 | 11.31 |
| `0510_pred_alone` | contrast | predicted NbNi2B2 alone | 38.36 |

## Database entries

Chemical system B-Nb-Ni. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/B-Nb-Ni.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/B2NbNi2_1378348.cif`
- search-match structures: `../../databases/searchmatch/<source>/B-Nb-Ni/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
