# YCoSi

Table 1 row 7, candidate list C. Discussed in the Supplement, Sec. S7, under YCoSi.

Outcome: Known Pnma YCoSi polymorph (+32 meV) formed, textured; not the on-hull P4/nmm structure

Category: Above-hull polymorph (probable). Structure previously reported at this composition.

P_disorder = 0.13. Database structures below the hull that are dynamically stable (MatterSim): 2 of 2.

## Buttons

Each folder in `buttons/` holds the synthesis record (`SYNTHESIS`, `SYNTHESIS.txt`, `STATUS`) and the measured scans.

| Button | Role | Mass loss (%) | Scans |
|---|---|---|---|
| 0503 | earlier attempt | 0.002 | `0503_HM_CoSiY_circular.txt` (Siemens D500, as cast) |
| 0538 | target | 0.039 | `0538_HM_CoSiY_circular.txt` (Siemens D500, as cast) |

Scan of record: `buttons/0538/0538_HM_CoSiY_circular.txt` (Siemens D500, powder, PMMA slide (inferred), as cast).

## Fits

Each folder in `fits/` is one GSAS-II refinement: `fit.txt` (observed and calculated pattern), `summary.json`, the GSAS-II project `fit.gpx` (in `work/` for most fits), the starting CIFs in `start/`, and `command.txt` to rerun it (see `../README.md`).

| Folder | Role | Model | R_wp (%) |
|---|---|---|---|
| `0538_FoR` | fit of record | Pnma YCoSi (known) + Y3Co2Si3 + YCo2Si2 | 10.09 |
| `0538_pred_P4nmm` | contrast | predicted P4/nmm in place of Pnma | 40.36 |

## Database entries

Chemical system Co-Si-Y. The entries are kept per system and shared with the other targets in it:

- hull energies: `../../databases/snapshot/<database>/Co-Si-Y.json` (oqmd_v1.7, mp, alexandria_pbe_2025.07.02, alexandria_pbesol_optimade), with the binary and elemental subsystems in the same folders
- predicted structure (the OQMD entry the target was chosen from): `../../databases/snapshot/structures/oqmd/CoSiY_1368530.cif`
- search-match structures: `../../databases/searchmatch/<source>/Co-Si-Y/` (oqmd_ternary_phases, mp_ternary_phases, alexandria_pbe_ternary_phases, cod)
