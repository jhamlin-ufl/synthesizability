# scripts/test_composition_deviation.py
"""
Test composition deviation parsing on all SYNTHESIS files.
Prints measured vs expected mole fractions and deviations.
"""

from pathlib import Path
import re
from synthesizability.parsers.synthesis import parse_synthesis_file

raw_dir = Path("data/raw")
synthesis_files = sorted(raw_dir.glob("*/SYNTHESIS"))

print(f"Testing {len(synthesis_files)} SYNTHESIS files\n")
print("=" * 70)

n_ok = 0
n_failed = 0

for sf in synthesis_files:
    folder = sf.parent.name
    parts = folder.split('_', 2)
    formula = parts[2] if len(parts) > 2 else None

    content = sf.read_text()
    result = parse_synthesis_file(content, formula)

    max_dev = result['composition_max_deviation']
    euc_dev = result['composition_euclidean_deviation']
    measured = result['composition_measured_fractions']
    expected = result['composition_expected_fractions']

    if max_dev is None:
        print(f"FAILED  {folder}")
        n_failed += 1
    else:
        flag = " *** HIGH ***" if max_dev > 0.02 else ""
        print(f"OK      {folder}  max={max_dev:.4f}  euc={euc_dev:.4f}{flag}")
        if measured and expected:
            for el in sorted(expected.keys()):
                exp = expected[el]
                meas = measured[el]
                diff = meas - exp
                print(f"         {el:2s}  expected={exp:.4f}  measured={meas:.4f}  diff={diff:+.4f}")
        n_ok += 1

print("=" * 70)
print(f"\nParsed OK: {n_ok}/{len(synthesis_files)}")
print(f"Failed:    {n_failed}/{len(synthesis_files)}")