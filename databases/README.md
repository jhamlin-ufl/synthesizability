# Databases

Entries of three public DFT databases and of the Crystallography Open Database for the 22
ternary chemical systems of the 34 targets and all their binary and elemental subsystems. The
analysis runs offline from these copies. Files are named by chemical system, for example
`Cu-Pd-Sc.json` (elements in alphabetical order). One system serves every target in it: Cu-Pd-Sc
serves ScPdCu4, ScPd2Cu, Sc2PdCu and Sc3Pd2Cu.

## `snapshot/`: hull energies, Table S2 and the phonon screen

| Path | Source | Retrieved |
|---|---|---|
| `oqmd_v1.7/<system>.json` | OQMD v1.7, local copy of the database download | 2026-03-02 |
| `mp/<system>.json` | Materials Project API (database version 2026.04.13 or earlier) | 2026-04-07 |
| `alexandria_pbe_2025.07.02/<system>.json` | Alexandria PBE release 2025.07.02 | 2026-10-01 |
| `alexandria_pbesol_optimade/<system>.json` | Alexandria PBEsol, OPTIMADE API | 2026-09-28 |
| `alexandria_pbesol_2023.12.29_list_c.csv` | the Alexandria PBEsol hull file that List C was selected from | 2026-09-28 |
| `structures/oqmd/` | the OQMD structure each target was selected from | 2026-02-23 |
| `structures/mp/`, `structures/alexandria_pbe/` | MP and Alexandria structures at the target compositions (index in `alexandria_pbe_index.csv`) | as above |
| `targets.csv` | the 34 targets with their chemical system and OQMD entry | |
| `MANIFEST.csv` | the source and SHA-256 of every file | |

Each `<system>.json` has `space`, `elements` and `entries`. Each entry has an id (`entry_id`, or
`mp_id`), `composition_id` (for example `"Al2 Au6 Mn1"`), `delta_e` (formation energy, eV/atom),
`stability` (the database's stored value, not used) and `icsd` (whether the entry was computed
from an ICSD record).

`workflow/scripts/hull_energies.py` recomputes each target's distance from the hull with
pymatgen and writes `results/hull_energies.csv`. The stored `stability` values are not used. For
Hf2MoIr the value OQMD stores (-18.6 meV/atom) disagrees with its own entries (-12.0 meV/atom).

MP holds one entry at ScPd2Au, mp-1096660, with a formation energy of +1.26 eV/atom. That entry
gives the MP value of +2078 meV/atom in the tables.

The PBEsol column of `results/hull_energies.csv` uses the OPTIMADE pull, which holds more entries
than the 2023.12.29 hull file. The two agree for seven of the eight List C targets. For YCoSi3 the
OPTIMADE hull gives -28.4 meV/atom and the hull file -44.3.

## `searchmatch/`: inputs of the automated search-match

| Path | Source |
|---|---|
| `oqmd_ternary_phases/<system>/`, `mp_ternary_phases/<system>/`, `alexandria_pbe_ternary_phases/<system>/` | the database entries with structures |
| `cod/<system>/` | Crystallography Open Database structures, with an index |
| `predicted_structures/<NNNN>_<formula>/` | the predicted structure of each target |
| `partially_disordered/<system>/` | ordered entries with some sites mixed, relaxed with MatterSim (our calculations) |
| `oqmd_entry_origins.tsv` | OQMD prototype of each entry, used to group structure types |

The code is in `workflow/src/synthesizability/phaseid/`.

## Licences

OQMD, the Materials Project and Alexandria are CC BY 4.0, and the COD is CC0. The extracts keep
these licences. Please cite:

- OQMD: J. E. Saal et al., JOM 65 (2013), doi:10.1007/s11837-013-0755-4; S. Kirklin et al., npj Comput. Mater. 1 (2015), doi:10.1038/npjcompumats.2015.10.
- Materials Project: A. Jain et al., APL Mater. 1 (2013), doi:10.1063/1.4812323.
- Alexandria: J. Schmidt et al., Mater. Today Phys. 48 (2024), doi:10.1016/j.mtphys.2024.101560; T. Cavignac et al., J. Phys. Mater. 9 (2026), doi:10.1088/2515-7639/ae6620.
- COD: S. Gražulis et al. (2012).
- List C used the PBEsol data of K. Gao et al., Nat. Commun. 16 (2025), doi:10.1038/s41467-025-63702-w. Only the hull entries are reproduced here.

No ICDD or ICSD files are included. Entries flagged `icsd` are the databases' own computed
entries.
