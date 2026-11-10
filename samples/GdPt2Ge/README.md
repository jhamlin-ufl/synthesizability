# GdPt2Ge

Table 1 row 3, candidate list B. Discussed in the Supplement, Sec. S7, under GdPt2Ge.

Outcome: YPd2Si-type Pnma (lowest in Alexandria), 94 wt% after anneal, small oxide fraction; OQMD Cmcm structure absent; new compound

Category: On-hull structure (probable). Not previously reported.

P_disorder = 0.18. Database structures below the hull that are dynamically stable (MatterSim): 2 of 2.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0480 | target | 0.33 | `0480_HM_GdGePt2.txt` (Siemens D500, as cast); `0480_HM_GdGePt2.xy` (PANalytical X'Pert MPD, as cast); `0480_HM_GdPt2Ge_annealed_ZB.txt` (Siemens D500, annealed 1000 C 1 week) |

Scan used: `buttons/0480/0480_HM_GdPt2Ge_annealed_ZB.txt` (Siemens D500, powder, zero-background Si plate, annealed 1000 C for one week).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0480_pnma_ox_genstrain` | main fit | GdPt2Ge annealed, Pnma + Gd2O3 | 11.63 |
| `0480_pnma` | contrast | GdPt2Ge annealed, Pnma alone, refined the same way as the Cmcm model | 14.17 |
| `0480_cmcm` | contrast | GdPt2Ge annealed, OQMD Cmcm alone | 41.48 |
| `0480_ascast_90321` | as-cast scan | GdPt2Ge as cast | 6.75 |
| `0480_structure` | quoted in the Supplement text, not plotted | main fit plus refined U_iso per site; Supplementary Tables of refined structures, with CIFs | 11.63 |

## Database entries

Chemical system Gd-Ge-Pt. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Gd-Ge-Pt.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/GdGePt2_1365099.cif`
- search-match structures: `../../databases/searchmatch/<source>/Gd-Ge-Pt/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
