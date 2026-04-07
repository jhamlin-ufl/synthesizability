# scripts/extract_mp_cifs.py
"""
Download CIF files from the Materials Project for all entries at the target
composition for each sample.

Reads per-space JSONs from data/external/mp_ternary_phases/<space>.json,
finds entries whose composition matches each sample's target formula, and
downloads their structures via the MP API.

CIF files are stored in:
  data/external/mp_ternary_phases/<space>/cifs/<formula>_<mp_id>_stab+<meV>meV.cif

Only entries at the exact target composition are downloaded (not all entries in
the space), keeping the download manageable.
"""
import json
import os
import re
import time
from itertools import combinations
from pathlib import Path
import sys

import pandas as pd
from mp_api.client import MPRester
from pymatgen.io.cif import CifWriter

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from synthesizability.oqmd import parse_elements_from_formula

DATA_DIR = Path("data/external/mp_ternary_phases")


def _get_api_key() -> str:
    key = os.environ.get("MP_API_KEY", "")
    if key:
        return key
    bashrc = Path.home() / ".bashrc"
    if bashrc.exists():
        for line in bashrc.read_text().splitlines():
            m = re.match(r"export\s+MP_API_KEY=(.+)", line.strip())
            if m:
                return m.group(1).strip().strip('"').strip("'")
    raise RuntimeError(
        "MP_API_KEY not found in environment or ~/.bashrc.\n"
        "Set it with: export MP_API_KEY=<your_key>"
    )


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
    if total == 0:
        return counts
    return {k: v / total for k, v in counts.items()}


def _fracs_match(entry_comp_id: str, target_comp_id: str, tol: float = 0.005) -> bool:
    """Check if two composition_id strings describe the same composition fractions."""
    ef = _to_fracs(_parse_composition(entry_comp_id))
    tf = _to_fracs(_parse_composition(target_comp_id))
    all_els = set(ef) | set(tf)
    return all(abs(ef.get(el, 0.0) - tf.get(el, 0.0)) < tol for el in all_els)


def _make_cif_filename(composition_id: str, mp_id: str, stability: float | None) -> str:
    """
    Build a CIF filename encoding composition, MP ID, and stability.
    e.g. 'Ta2Ti1Mo1_mp-530936_stab+0meV.cif'
    """
    compact = composition_id.replace(" ", "")
    if stability is None:
        stab_str = "stabNone"
    else:
        meV = round(stability * 1000)
        stab_str = f"stab+{meV}meV" if meV >= 0 else f"stab{meV}meV"
    safe_mp_id = mp_id.replace("/", "_")
    return f"{compact}_{safe_mp_id}_{stab_str}.cif"


def get_target_entries(formula: str) -> list[tuple[str, dict]]:
    """
    Find all MP entries whose composition matches the sample formula.
    Returns list of (space, entry) tuples.
    """
    from synthesizability.oqmd import parse_formula_to_oqmd
    target_comp_id = parse_formula_to_oqmd(formula)  # e.g. 'Mo1 Ta2 Ti1'

    elements = parse_elements_from_formula(formula)
    matches = []

    for r in range(1, len(elements) + 1):
        for combo in combinations(elements, r):
            space = "-".join(sorted(combo))
            json_path = DATA_DIR / f"{space}.json"
            if not json_path.exists():
                continue
            payload = json.loads(json_path.read_text())
            for entry in payload["entries"]:
                if _fracs_match(entry["composition_id"], target_comp_id):
                    matches.append((space, entry))

    return matches


def main():
    api_key = _get_api_key()
    print(f"MP API key: {api_key[:8]}...")

    print("Loading synthesis data...")
    df = pd.read_csv("data/processed/synthesis_data_no_disorder.csv")
    formulas = df["formula"].dropna().unique().tolist()
    print(f"Found {len(formulas)} unique formulas\n")

    with MPRester(api_key) as mpr:
        version = mpr.get_database_version()
        print(f"Connected to MP database version: {version}\n")

        grand_downloaded = 0
        grand_skipped = 0
        grand_failed = 0

        for formula in sorted(formulas):
            matches = get_target_entries(formula)
            if not matches:
                print(f"{formula}: no MP entries at target composition")
                continue

            print(f"{formula}: {len(matches)} target-composition entries")

            for space, entry in matches:
                mp_id = entry["mp_id"]
                comp_id = entry["composition_id"]
                stability = entry["stability"]

                cif_dir = DATA_DIR / space / "cifs"
                cif_dir.mkdir(parents=True, exist_ok=True)

                filename = _make_cif_filename(comp_id, mp_id, stability)
                cif_path = cif_dir / filename

                if cif_path.exists():
                    print(f"  {mp_id} ({comp_id}): skipped (exists)")
                    grand_skipped += 1
                    continue

                try:
                    structure = mpr.get_structure_by_material_id(mp_id)
                    writer = CifWriter(structure)
                    writer.write_file(str(cif_path))
                    stab_meV = round(stability * 1000) if stability is not None else None
                    print(f"  {mp_id} ({comp_id}): downloaded (stab={stab_meV} meV/atom)")
                    grand_downloaded += 1
                    time.sleep(0.1)  # be polite to the API
                except Exception as e:
                    print(f"  {mp_id} ({comp_id}): FAILED — {e}")
                    grand_failed += 1

    print(f"\n{'='*60}")
    print(f"CIF extraction complete")
    print(f"  Downloaded: {grand_downloaded}")
    print(f"  Skipped:    {grand_skipped} (already existed)")
    print(f"  Failed:     {grand_failed}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
