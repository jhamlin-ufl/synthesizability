# Ti2Pt8Ge

Table 1 row 33, candidate list B. Discussed in the Supplement, Sec. S7, under Ti2Pt8Ge.

Outcome: Pt-rich fcc-type phase plus unidentified lines; on-hull structure absent

Category: On-hull structure did not form but structure(s) unresolved (probable). Not previously reported.

P_disorder = 0.26. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0492 | target | 0.53 | `0492_HM_Ti2GePt8_circular.txt` (Siemens D500, as cast) |

Scan used: `buttons/0492/0492_HM_Ti2GePt8_circular.txt` (Siemens D500, bulk chunk, as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0492_best` | best fit found (not a satisfactory description; not shown in the Supplement) | Pt-rich fcc-derived Ti3Pt7 type (bulk chunk) | 14.92 |
| `0492_pred` | contrast | ordered on-hull P4/mmm alone | 22.47 |

## Database entries

Chemical system Ge-Pt-Ti. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Ge-Pt-Ti.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/Ti2GePt8_2149411.cif`
- search-match structures: `../../databases/searchmatch/<source>/Ge-Pt-Ti/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
