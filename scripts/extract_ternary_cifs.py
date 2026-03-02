# scripts/extract_ternary_cifs.py
"""
Extract CIF files for all phases in each chemical space JSON.
Reads JSONs from data/external/oqmd_ternary_phases/<space>.json and writes
CIFs to data/external/oqmd_ternary_phases/<space>/cifs/.

CIF filenames encode composition, entry ID, and stability in meV/atom.
Example: Co2Si2Y1_647098_stab-92meV.cif
"""
import json
import time
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from synthesizability.oqmd import (
    check_database_exists,
    get_structure_from_db,
    make_cif_filename,
)
from pymatgen.io.cif import CifWriter

DATA_DIR = Path("data/external/oqmd_ternary_phases")


def process_space(json_path: Path) -> tuple[int, int, int]:
    """
    Extract CIFs for one chemical space.

    Returns:
        (n_extracted, n_skipped, n_failed) counts
    """
    payload = json.loads(json_path.read_text())
    entries = payload["entries"]
    space = payload["space"]

    if not entries:
        return 0, 0, 0

    cif_dir = DATA_DIR / space / "cifs"
    cif_dir.mkdir(parents=True, exist_ok=True)

    n_extracted = 0
    n_skipped = 0
    n_failed = 0

    for entry in entries:
        entry_id = entry["entry_id"]
        composition_id = entry["composition_id"]
        stability = entry["stability"]

        filename = make_cif_filename(composition_id, entry_id, stability)
        cif_path = cif_dir / filename

        if cif_path.exists():
            n_skipped += 1
            continue

        structure = get_structure_from_db(entry_id)

        if structure is None:
            print(f"    ⚠ Failed: {composition_id} (entry {entry_id})")
            n_failed += 1
            continue

        writer = CifWriter(structure)
        writer.write_file(str(cif_path))
        n_extracted += 1

    return n_extracted, n_skipped, n_failed


def main():
    if not check_database_exists():
        print("ERROR: OQMD database not found!")
        print("Please run: poetry run python scripts/validate_oqmd_database.py")
        sys.exit(1)

    json_paths = sorted(DATA_DIR.glob("*.json"))
    if not json_paths:
        print(f"No JSON files found in {DATA_DIR}")
        print("Run scripts/query_ternary_phases.py first.")
        sys.exit(1)

    # Skip empty spaces upfront
    nonempty = [p for p in json_paths if json.loads(p.read_text())["entries"]]
    print(f"Found {len(json_paths)} space JSONs, {len(nonempty)} non-empty")
    print(f"Output: {DATA_DIR}/<space>/cifs/\n")

    total_start = time.time()
    grand_extracted = 0
    grand_skipped = 0
    grand_failed = 0

    for i, json_path in enumerate(nonempty):
        space = json_path.stem
        t0 = time.time()
        print(f"[{i+1}/{len(nonempty)}] {space}...", end=" ", flush=True)

        n_extracted, n_skipped, n_failed = process_space(json_path)
        elapsed = time.time() - t0

        print(f"extracted={n_extracted} skipped={n_skipped} failed={n_failed} ({elapsed:.1f}s)")

        grand_extracted += n_extracted
        grand_skipped += n_skipped
        grand_failed += n_failed

    total_elapsed = time.time() - total_start
    print(f"\n{'='*60}")
    print(f"Extraction complete in {total_elapsed:.1f}s")
    print(f"  Extracted: {grand_extracted}")
    print(f"  Skipped:   {grand_skipped} (already existed)")
    print(f"  Failed:    {grand_failed}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()