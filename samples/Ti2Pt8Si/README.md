# Ti2Pt8Si

Table 1 row 23, candidate list B. Discussed in the Supplement, Sec. S7, under Ti2Pt8Si.

Outcome: TiPt5Si-type and Pt-rich fcc phase formed; on-hull structure absent

Category: Multiphase, on-hull structure not formed (probable). Not previously reported.

P_disorder = 0.19. Database structures below the hull that are dynamically stable (MatterSim): 1 of 1.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0494 | target | 0.44 | `0494_HM_Ti2SiPt8.txt` (Siemens D500, as cast) |

Scan of record: `buttons/0494/0494_HM_Ti2SiPt8.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0494_FoR_hs` | fit of record | TiSiPt5-type + fcc-derived Ti-Pt | 7.86 |
| `0494_pred_alone` | contrast | predicted Ti2Pt8Si alone | 25.18 |
| `0494_pred_for_TiSiPt5_hs` | quoted in the Supplement text, not plotted | predicted Ti2Pt8Si in place of the TiPt5Si-type phase | 14.50 |

## Database entries

Chemical system Pt-Si-Ti. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Pt-Si-Ti.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/Ti2SiPt8_2149228.cif`
- search-match structures: `../../databases/searchmatch/<source>/Pt-Si-Ti/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
