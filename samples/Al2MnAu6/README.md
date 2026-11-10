# Al2MnAu6

Table 1 row 15, candidate list B. Discussed in the Supplement, Sec. S7, under Al2MnAu6.

Outcome: C11b Au2(Al,Mn), partial Au order; Mn/Al order undetermined

Category: Disordered parent of on-hull structure (probable). Not previously reported.

P_disorder = 0.62. Database structures below the hull that are dynamically stable (MatterSim): 0 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0479 | target | 2.01 | `0479_HM_MnAl2Au6.txt` (Siemens D500, as cast); `0479_HM_MnAl2Au6.xy` (PANalytical X'Pert MPD, as cast) |

Scan used: `buttons/0479/0479_HM_MnAl2Au6.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0479_FoR_C11b_pos_hs` | main fit | C11b Au2(Al,Mn), Au order eta 0.5, Au z free | 6.18 |
| `0479_pred_C2m_eta0.50_pos_hs` | contrast | C2/m with predicted Mn/Al order, same Au disorder, positions free | 5.93 |

## Database entries

Chemical system Al-Au-Mn. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Al-Au-Mn.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/MnAl2Au6_1992237.cif`
- search-match structures: `../../databases/searchmatch/<source>/Al-Au-Mn/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
