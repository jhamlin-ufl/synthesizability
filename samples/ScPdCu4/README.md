# ScPdCu4

Table 1 row 34, candidate list A. Discussed in the Supplement, Sec. S7, under ScPdCu4.

Outcome: Broad pattern, no well-defined lines; on-hull structure absent

Category: On-hull structure did not form but structure(s) unresolved (tentative). Not previously reported.

P_disorder = 0.25. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0472 | target | 5.06 | `0472_HM_ScCu4Pd.txt` (Siemens D500, as cast); `0472_HM_ScCu4Pd.xy` (PANalytical X'Pert MPD, as cast) |
| 0558 | remake | 13.04 | none |
| 0566 | remake | 4.72 | `0566_HM_ScCu4Pd.txt` (Siemens D500, as cast) |
| 0613 | remake | 4.48 | `0613_HM_ScCu4Pd.txt` (Siemens D500, as cast) |

Scan of record: `buttons/0472/0472_HM_ScCu4Pd.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0472_FoR_hs` | best fit found (not a satisfactory description; not shown in the Supplement) | Cu fcc + B2 (Cu,Pd)(Sc) + ScPd | 5.94 |
| `0472_pred_alone` | contrast | predicted ScPdCu4 alone | 8.10 |

## Database entries

Chemical system Cu-Pd-Sc. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Cu-Pd-Sc.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/ScCu4Pd_1991578.cif`
- search-match structures: `../../databases/searchmatch/<source>/Cu-Pd-Sc/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
