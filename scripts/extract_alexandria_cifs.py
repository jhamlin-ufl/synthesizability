# scripts/extract_alexandria_cifs.py
"""
Fetch CIF files from Alexandria for all entries at the target composition
for each sample.  Handles both PBE and PBEsol in one run.

Reads per-space JSONs from:
  data/external/alexandria_pbe_ternary_phases/<space>.json
  data/external/alexandria_pbesol_ternary_phases/<space>.json

Finds entries whose composition matches each sample's target formula, then
fetches full structures via the Alexandria OPTIMADE single-entry endpoint
and writes CIF files using pymatgen.

CIF files are stored in:
  data/external/alexandria_pbe_ternary_phases/<space>/cifs/
      <compact_comp>_<entry_id>_stab+<meV>meV.cif
  data/external/alexandria_pbesol_ternary_phases/<space>/cifs/
      (same layout)
"""
import json
import re
import time
from itertools import combinations
from pathlib import Path
import sys

import pandas as pd
import requests
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from synthesizability.oqmd import parse_elements_from_formula, parse_formula_to_oqmd

FUNCTIONALS = {
    "pbe":    ("data/external/alexandria_pbe_ternary_phases",
               "https://alexandria.icams.rub.de/pbe/v1"),
    "pbesol": ("data/external/alexandria_pbesol_ternary_phases",
               "https://alexandria.icams.rub.de/pbesol/v1"),
}

STRUCTURE_FIELDS = (
    "lattice_vectors,cartesian_site_positions,"
    "species,species_at_sites,nsites"
)
SLEEP_BETWEEN = 0.5


def _parse_composition(composition_id: str) -> dict[str, float]:
    """'Mo2 Ta1 Ti2' → {'Mo': 2.0, 'Ta': 1.0, 'Ti': 2.0}"""
    out = {}
    for token in composition_id.split():
        m = re.match(r"([A-Z][a-z]*)(\d+(?:\.\d+)?)", token)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def _to_fracs(counts: dict[str, float]) -> dict[str, float]:
    total = sum(counts.values())
    return {k: v / total for k, v in counts.items()} if total else counts


def _fracs_match(entry_comp_id: str, target_comp_id: str, tol: float = 0.005) -> bool:
    ef = _to_fracs(_parse_composition(entry_comp_id))
    tf = _to_fracs(_parse_composition(target_comp_id))
    all_els = set(ef) | set(tf)
    return all(abs(ef.get(el, 0.0) - tf.get(el, 0.0)) < tol for el in all_els)


def _make_cif_filename(composition_id: str, entry_id: str,
                        stability: float | None) -> str:
    compact = composition_id.replace(" ", "")
    if stability is None:
        stab_str = "stabNone"
    else:
        meV = round(stability * 1000)
        stab_str = f"stab+{meV}meV" if meV >= 0 else f"stab{meV}meV"
    return f"{compact}_{entry_id}_{stab_str}.cif"


def _fetch_structure(base_url: str, entry_id: str) -> Structure | None:
    """Fetch a single structure from the OPTIMADE single-entry endpoint."""
    url = f"{base_url}/structures/{entry_id}"
    params = {"response_fields": STRUCTURE_FIELDS}
    try:
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        print(f"    ERROR fetching {entry_id}: {e}")
        return None

    attrs = data.get("data", {}).get("attributes", {})

    lattice_vectors = attrs.get("lattice_vectors")
    positions       = attrs.get("cartesian_site_positions")
    species_defs    = attrs.get("species", [])
    species_at_sites = attrs.get("species_at_sites", [])

    if not lattice_vectors or not positions or not species_at_sites:
        print(f"    WARNING: incomplete structure data for {entry_id}")
        return None

    # Build species name → chemical symbol map
    # species_defs is a list of {"name": "Al", "chemical_symbols": ["Al"], ...}
    species_map = {}
    for sp in species_defs:
        name = sp.get("name", "")
        symbols = sp.get("chemical_symbols", [])
        # Take the first non-vacancy symbol
        for sym in symbols:
            if sym != "X":
                species_map[name] = sym
                break

    try:
        lattice   = Lattice(lattice_vectors)
        site_species = [species_map.get(s, s) for s in species_at_sites]
        structure = Structure(
            lattice, site_species, positions, coords_are_cartesian=True
        )
        return structure
    except Exception as e:
        print(f"    ERROR building Structure for {entry_id}: {e}")
        return None


def get_target_entries(formula: str, data_dir: Path) -> list[tuple[str, dict]]:
    """Find all entries whose composition matches the sample formula."""
    target_comp_id = parse_formula_to_oqmd(formula)
    elements = parse_elements_from_formula(formula)
    matches = []

    for r in range(1, len(elements) + 1):
        for combo in combinations(elements, r):
            space = "-".join(sorted(combo))
            json_path = data_dir / f"{space}.json"
            if not json_path.exists():
                continue
            payload = json.loads(json_path.read_text())
            for entry in payload["entries"]:
                if _fracs_match(entry["composition_id"], target_comp_id):
                    matches.append((space, entry))

    return matches


def main():
    print("Loading synthesis data...")
    df = pd.read_csv("data/processed/synthesis_data_no_disorder.csv")
    formulas = df["formula"].dropna().unique().tolist()
    print(f"Found {len(formulas)} unique formulas\n")

    for functional, (data_dir_str, base_url) in FUNCTIONALS.items():
        data_dir = Path(data_dir_str)

        print(f"{'='*60}")
        print(f"Extracting CIFs — Alexandria {functional.upper()}")
        print(f"Data dir: {data_dir}")
        print(f"{'='*60}")

        grand_downloaded = grand_skipped = grand_failed = grand_nomatch = 0

        for formula in sorted(formulas):
            matches = get_target_entries(formula, data_dir)
            if not matches:
                print(f"{formula}: no Alexandria {functional.upper()} entries at target composition")
                grand_nomatch += 1
                continue

            print(f"{formula}: {len(matches)} target-composition entries")

            for space, entry in matches:
                entry_id  = entry["entry_id"]
                comp_id   = entry["composition_id"]
                stability = entry["stability"]

                cif_dir = data_dir / space / "cifs"
                cif_dir.mkdir(parents=True, exist_ok=True)

                filename = _make_cif_filename(comp_id, entry_id, stability)
                cif_path = cif_dir / filename

                if cif_path.exists():
                    print(f"  {entry_id} ({comp_id}): skipped (exists)")
                    grand_skipped += 1
                    continue

                structure = _fetch_structure(base_url, entry_id)
                if structure is None:
                    grand_failed += 1
                    continue

                try:
                    CifWriter(structure).write_file(str(cif_path))
                    stab_meV = round(stability * 1000) if stability is not None else None
                    print(f"  {entry_id} ({comp_id}): downloaded (stab={stab_meV} meV/atom)")
                    grand_downloaded += 1
                except Exception as e:
                    print(f"  {entry_id} ({comp_id}): FAILED writing CIF — {e}")
                    grand_failed += 1

                time.sleep(SLEEP_BETWEEN)

        print(f"\nDone with {functional.upper()}")
        print(f"  Downloaded: {grand_downloaded}")
        print(f"  Skipped:    {grand_skipped} (already existed)")
        print(f"  No match:   {grand_nomatch}")
        print(f"  Failed:     {grand_failed}")
        print()


if __name__ == "__main__":
    main()
