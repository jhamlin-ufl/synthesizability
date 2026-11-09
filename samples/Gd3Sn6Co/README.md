# Gd3Sn6Co

Table 1 row 18, candidate list B. Discussed in the Supplement, Sec. S7, under Gd3Sn6Co.

Outcome: Known CeNiSi2-type GdCo_xSn2; Co ordering of the on-hull structure absent

Category: Disordered parent of on-hull structure (probable). Structure previously reported at this composition.

P_disorder = 0.24. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0489 | earlier attempt | 4.54 | `0489_HM_Gd3CoSn6_circular.txt` (Siemens D500, as cast) |
| 0536 | remake | 9.46 | `0536_HM_Gd3CoSn6_circular.txt` (Siemens D500, as cast) |
| 0556 | target | 2.04 | `0556_HM_Gd3CoSn6.txt` (Siemens D500, as cast) |

Scan of record: `buttons/0556/0556_HM_Gd3CoSn6.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0556_FoR` | fit of record | CeNiSi2-type GdCo_xSn_2 (Gil 2004), partners below 1 wt% | 13.39 |
| `0556_pred_c_hs` | contrast | predicted Co-ordered Cmcm Gd3CoSn6 in place of GdCo_xSn_2, sample displacement held | 13.26 |

## Database entries

Chemical system Co-Gd-Sn. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Co-Gd-Sn.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/Gd3CoSn6_1790943.cif`
- search-match structures: `../../databases/searchmatch/<source>/Co-Gd-Sn/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
