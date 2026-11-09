# TiPt5Ge

Table 1 row 22, candidate list B. Discussed in the Supplement, Sec. S7, under TiPt5Ge.

Outcome: Pt + Fe2P-type Pt2Ge; on-hull structure absent

Category: Multiphase, on-hull structure not formed (probable). Not previously reported.

P_disorder = 0.20. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0491 | target | 1.40 | `0491_HM_TiGePt5_circular.txt` (Siemens D500, as cast) |

Scan of record: `buttons/0491/0491_HM_TiGePt5_circular.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0491_FoR` | fit of record | Pt + Pt2Ge | 5.36 |
| `0491_pred_alone` | contrast | predicted TiPt5Ge alone | 22.72 |

## Database entries

Chemical system Ge-Pt-Ti. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Ge-Pt-Ti.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/TiGePt5_2146080.cif`
- search-match structures: `../../databases/searchmatch/<source>/Ge-Pt-Ti/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
