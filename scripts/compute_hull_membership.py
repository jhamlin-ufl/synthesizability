# scripts/compute_hull_membership.py
"""
For each unique sample formula, determine which phase databases have
that composition on (or below) the convex hull.

Outputs data/processed/hull_membership.csv with columns:
    formula            - sample formula string
    hull_sources       - Python list repr of database names on hull
    hull_source_count  - integer count (0–4)
"""
import sys
import json
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from synthesizability.dashboard_plugins.ternary_phases import get_hull_sources


def main():
    formulas_path = Path("data/processed/synthesis_data_no_disorder.csv")
    output_path   = Path("data/processed/hull_membership.csv")

    df_in = pd.read_csv(formulas_path)
    formulas = sorted(df_in["formula"].dropna().unique())
    print(f"Computing hull membership for {len(formulas)} unique formulas...")

    source_names = ["OQMD", "MP", "Alexandria PBE", "Alexandria PBEsol"]
    col_names    = ["hull_oqmd", "hull_mp", "hull_alex_pbe", "hull_alex_pbesol"]

    rows = []
    for formula in formulas:
        sources = get_hull_sources(formula)
        row = {
            "formula":           formula,
            "hull_sources":      str(sources),   # Python list repr, consistent with other list cols
            "hull_source_count": len(sources),
        }
        for src, col in zip(source_names, col_names):
            row[col] = src in sources
        rows.append(row)
        if sources:
            print(f"  {formula}: on hull in {sources}")

    pd.DataFrame(rows).to_csv(output_path, index=False)
    n_on = sum(1 for r in rows if r["hull_source_count"] > 0)
    print(f"\n✓ Wrote {len(rows)} rows to {output_path}")
    print(f"  {n_on}/{len(rows)} formulas on hull in at least one database")


if __name__ == "__main__":
    main()
