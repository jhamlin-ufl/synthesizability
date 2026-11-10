# PtNiGe2

Table 1 row 17, candidate list B. Discussed in the Supplement, Sec. S7, under PtNiGe2.

Outcome: MnP-type (Ni,Pt)Ge, parent of the on-hull structure, Ni/Pt order absent; pattern unresolved

Category: Disordered parent of on-hull structure (tentative). Not previously reported.

P_disorder = 0.51. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0481 | target | 0.25 | `0481_HM_NiGe2Pt.txt` (Siemens D500, as cast); `0481_HM_NiGe2Pt.xy` (PANalytical X'Pert MPD, as cast) |

Scan used: `buttons/0481/0481_HM_NiGe2Pt.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0481_FoR_w14` | main fit | MnP-type (Ni0.5Pt0.5)Ge, Ni/Pt mixed, alone | 21.25 |
| `0481_pred_w14` | contrast | ordered P2_1/m prediction alone | 36.31 |

## Database entries

Chemical system Ge-Ni-Pt. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Ge-Ni-Pt.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/NiGe2Pt_1646030.cif`
- search-match structures: `../../databases/searchmatch/<source>/Ge-Ni-Pt/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
