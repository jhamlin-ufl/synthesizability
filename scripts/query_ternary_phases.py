# scripts/query_ternary_phases.py
"""
Query OQMD for all phases in each unique chemical space (unary, binary, ternary)
across all samples. Each space is queried once and stored as a JSON file at
data/external/oqmd_ternary_phases/<space>.json, e.g. Al-Co-Fe.json.

Spaces are derived from sample formulas in synthesis_data_no_disorder.csv.
"""
import json
import time
import pandas as pd
from itertools import combinations
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from synthesizability.oqmd import (
    check_database_exists,
    parse_elements_from_formula,
    query_exact_space,
)

OUTPUT_DIR = Path("data/external/oqmd_ternary_phases")


def get_all_spaces(formulas: list[str]) -> list[list[str]]:
    """
    Get all unique chemical spaces (as sorted element lists) across all formulas.
    Includes all subspaces: unaries, binaries, and ternaries.
    """
    space_set = set()
    for formula in formulas:
        elements = parse_elements_from_formula(formula)
        for r in range(1, len(elements) + 1):
            for combo in combinations(elements, r):
                space_set.add(tuple(sorted(combo)))
    return [list(s) for s in sorted(space_set)]


def space_to_key(elements: list[str]) -> str:
    """Convert element list to hyphen-separated string key, e.g. ['Al','Co','Fe'] -> 'Al-Co-Fe'."""
    return '-'.join(elements)


def main():
    if not check_database_exists():
        print("ERROR: OQMD database not found!")
        print("Please run: poetry run python scripts/validate_oqmd_database.py")
        sys.exit(1)

    print("Loading synthesis data...")
    df = pd.read_csv("data/processed/synthesis_data_no_disorder.csv")
    print(f"Found {len(df)} samples")

    all_spaces = get_all_spaces(df["formula"].tolist())
    print(f"Found {len(all_spaces)} unique chemical spaces")
    by_order = {n: sum(1 for s in all_spaces if len(s) == n) for n in [1, 2, 3]}
    for n, count in by_order.items():
        print(f"  {n}-element spaces: {count}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    total_start = time.time()
    n_total = len(all_spaces)
    n_success = 0
    n_empty = 0
    n_skipped = 0

    for idx, elements in enumerate(all_spaces):
        key = space_to_key(elements)
        out_path = OUTPUT_DIR / f"{key}.json"

        if out_path.exists():
            print(f"[{idx+1}/{n_total}] {key}: skipping (already exists)")
            n_skipped += 1
            continue

        print(f"[{idx+1}/{n_total}] {key}...", end=" ", flush=True)

        t0 = time.time()
        entries = query_exact_space(elements)
        elapsed = time.time() - t0

        if not entries:
            print(f"0 entries ({elapsed:.1f}s)")
            n_empty += 1
        else:
            print(f"{len(entries)} entries ({elapsed:.1f}s)")
            n_success += 1

        # Write JSON regardless (empty list is valid — avoids re-querying)
        payload = {
            "space": key,
            "elements": elements,
            "entries": entries,
        }
        out_path.write_text(json.dumps(payload, indent=2))

    total_elapsed = time.time() - total_start
    print(f"\n{'='*60}")
    print(f"Query complete in {total_elapsed:.1f}s")
    print(f"  With entries: {n_success}")
    print(f"  Empty:        {n_empty}")
    print(f"  Skipped:      {n_skipped} (already existed)")
    print(f"  Output:       {OUTPUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()