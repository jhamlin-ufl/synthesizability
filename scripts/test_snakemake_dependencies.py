#!/usr/bin/env python3
"""
Test Snakemake dependency triggers.

Checks that modifying files in various locations causes the expected rules
to be marked for rerun, and that unrelated rules are NOT triggered.

Uses --dry-run to detect planned reruns without executing anything.
Makes real content changes to trigger checksum-based rerun detection,
and restores files afterward.
"""

import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent


def run_dry() -> set[str]:
    """
    Run snakemake dry-run and return set of rule names that would execute.
    Parses the 'Job stats' table from dry-run output.
    """
    result = subprocess.run(
        ["poetry", "run", "snakemake", "--dry-run", "--cores", "1"],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    output = result.stdout + result.stderr

    rules = set()
    in_stats = False
    for line in output.splitlines():
        if line.strip().startswith("job ") and "count" in line:
            in_stats = True
            continue
        if in_stats:
            if line.strip().startswith("---"):
                continue
            if line.strip().startswith("total"):
                break
            parts = line.strip().split()
            if parts:
                rules.add(parts[0])

    rules.discard("all")
    return rules


def modify_file(path: Path) -> bytes:
    """Append a harmless comment to force a real checksum change. Returns original content."""
    original = path.read_bytes()
    path.write_bytes(original + b"\n# _test_trigger\n")
    return original


def restore_file(path: Path, original: bytes) -> None:
    """Restore file to original content."""
    path.write_bytes(original)


def check_baseline() -> bool:
    """Confirm nothing needs to run before we start."""
    rules = run_dry()
    if rules:
        print(f"  WARNING: baseline is not clean — {rules} would run")
        return False
    print("  OK: baseline is clean")
    return True


def run_test(
    description: str,
    modify_files: list[Path],
    expect_triggered: set[str],
    expect_clean: set[str],
) -> bool:
    """
    Modify files, check dry-run output, restore files, report result.
    Returns True if test passed.
    """
    print(f"\n{'─'*60}")
    print(f"TEST: {description}")
    print(f"  Modifying: {[str(f.relative_to(REPO_ROOT)) for f in modify_files]}")

    originals = {}
    for f in modify_files:
        originals[f] = modify_file(f)

    try:
        triggered = run_dry()
        print(f"  Triggered: {triggered}")

        missing = expect_triggered - triggered
        unexpected = expect_clean & triggered

        passed = not missing and not unexpected

        if missing:
            print(f"  FAIL: expected but not triggered: {missing}")
        if unexpected:
            print(f"  FAIL: triggered but should not have been: {unexpected}")
        if passed:
            print(f"  PASS")

        return passed

    finally:
        for f, content in originals.items():
            restore_file(f, content)


def main():
    print("=" * 60)
    print("Snakemake Dependency Trigger Tests")
    print("=" * 60)

    if not check_baseline():
        print("Baseline not clean — please run snakemake first to bring outputs up to date.")
        return

    results = []

    results.append(run_test(
        description="Modify dashboard plugin → only generate_dashboard",
        modify_files=[
            REPO_ROOT / "src/synthesizability/dashboard_plugins/supercon.py"
        ],
        expect_triggered={"generate_dashboard"},
        expect_clean={"compute_disorder_cache", "analyze_susceptibility",
                      "build_dataframe", "compute_supercon_cache"},
    ))

    results.append(run_test(
        description="Modify disorder.py → only compute_disorder_cache (and downstream)",
        modify_files=[
            REPO_ROOT / "src/synthesizability/disorder.py"
        ],
        expect_triggered={"compute_disorder_cache"},
        expect_clean={"analyze_susceptibility", "compute_supercon_cache"},
    ))

    results.append(run_test(
        description="Modify susceptibility.py → only analyze_susceptibility (and downstream)",
        modify_files=[
            REPO_ROOT / "src/synthesizability/susceptibility.py"
        ],
        expect_triggered={"analyze_susceptibility"},
        expect_clean={"compute_disorder_cache", "compute_supercon_cache",
                      "build_dataframe_for_formulas"},
    ))

    results.append(run_test(
        description="Modify oqmd.py → only oqmd rules",
        modify_files=[
            REPO_ROOT / "src/synthesizability/oqmd.py"
        ],
        expect_triggered={"query_oqmd_hulls", "extract_oqmd_structures"},
        expect_clean={"compute_disorder_cache", "analyze_susceptibility",
                      "compute_supercon_cache"},
    ))

    io_files = list((REPO_ROOT / "src/synthesizability/io").rglob("*.py"))
    if io_files:
        results.append(run_test(
            description="Modify io source → build_dataframe rules + analyze_susceptibility",
            modify_files=[io_files[0]],
            expect_triggered={"build_dataframe_for_formulas", "build_dataframe",
                              "analyze_susceptibility"},
            expect_clean=set(),
        ))

    chi_files = list((REPO_ROOT / "data/raw").rglob("*chiAC*.txt"))
    if chi_files:
        results.append(run_test(
            description="Modify chi data file → analyze_susceptibility + generate_dashboard",
            modify_files=[chi_files[0]],
            expect_triggered={"analyze_susceptibility", "generate_dashboard"},
            expect_clean={"compute_disorder_cache", "compute_supercon_cache",
                          "build_dataframe_for_formulas"},
        ))

    status_files = list((REPO_ROOT / "data/raw").rglob("STATUS"))
    if status_files:
        results.append(run_test(
            description="Modify STATUS file → build_dataframe_for_formulas + build_dataframe",
            modify_files=[status_files[0]],
            expect_triggered={"build_dataframe_for_formulas", "build_dataframe"},
            expect_clean={"analyze_susceptibility"},
        ))

    print(f"\n{'='*60}")
    passed = sum(results)
    total = len(results)
    print(f"Results: {passed}/{total} tests passed")
    if passed < total:
        print("SOME TESTS FAILED — review Snakefile dependencies")
    else:
        print("All tests passed ✓")
    print("=" * 60)


if __name__ == "__main__":
    main()
