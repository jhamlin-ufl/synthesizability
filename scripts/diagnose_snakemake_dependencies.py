# scripts/diagnose_snakemake_dependencies.py
"""
Diagnostic script to verify Snakemake file globs are capturing the right files.
"""

from pathlib import Path

def main():
    print("="*80)
    print("SNAKEMAKE DEPENDENCY DIAGNOSTIC")
    print("="*80)
    
    # Check raw data directory exists
    raw_dir = Path("data/raw")
    if not raw_dir.exists():
        print(f"\nERROR: {raw_dir} does not exist!")
        return
    
    # Get all files
    all_files = [p for p in raw_dir.rglob("*") if p.is_file()]
    print(f"\nTotal files in data/raw: {len(all_files)}")
    
    # Chi files
    chi_files = [p for p in raw_dir.rglob("*chiAC*.txt")]
    print(f"\nChi data files (*chiAC*.txt): {len(chi_files)}")
    if chi_files:
        print("  Examples:")
        for f in sorted(chi_files)[:3]:
            print(f"    {f.relative_to(raw_dir)}")
    
    # STATUS files
    status_files = [p for p in raw_dir.rglob("STATUS")]
    print(f"\nSTATUS files: {len(status_files)}")
    if status_files:
        print("  Examples:")
        for f in sorted(status_files)[:3]:
            print(f"    {f.relative_to(raw_dir)}")
    
    # SYNTHESIS files
    synthesis_files = [p for p in raw_dir.rglob("SYNTHESIS")]
    print(f"\nSYNTHESIS files: {len(synthesis_files)}")
    if synthesis_files:
        print("  Examples:")
        for f in sorted(synthesis_files)[:3]:
            print(f"    {f.relative_to(raw_dir)}")
    
    # XRD files (.xy)
    xrd_xy_files = [p for p in raw_dir.rglob("*.xy")]
    print(f"\nXRD files (*.xy): {len(xrd_xy_files)}")
    if xrd_xy_files:
        print("  Examples:")
        for f in sorted(xrd_xy_files)[:3]:
            print(f"    {f.relative_to(raw_dir)}")
    
    # Other .txt files (excluding chi)
    other_txt = [p for p in raw_dir.rglob("*.txt") if "chiAC" not in p.name]
    print(f"\nOther .txt files (excluding chiAC): {len(other_txt)}")
    if other_txt:
        print("  Examples:")
        for f in sorted(other_txt)[:5]:
            print(f"    {f.relative_to(raw_dir)}")
    
    # Check for unexpected file types
    all_extensions = set(p.suffix for p in all_files if p.suffix)
    print(f"\nAll file extensions found: {sorted(all_extensions)}")
    
    # Files without extensions
    no_ext = [p for p in all_files if not p.suffix]
    print(f"\nFiles without extension: {len(no_ext)}")
    if no_ext:
        unique_names = set(p.name for p in no_ext)
        print(f"  Unique names: {sorted(unique_names)}")
    
    # Coverage check
    dataframe_inputs = set(status_files + synthesis_files + xrd_xy_files + other_txt)
    chi_only = set(chi_files)
    covered = dataframe_inputs | chi_only
    uncovered = set(all_files) - covered
    
    print("\n" + "="*80)
    print("COVERAGE ANALYSIS")
    print("="*80)
    print(f"Files tracked for dataframe building: {len(dataframe_inputs)}")
    print(f"Files tracked for chi analysis: {len(chi_only)}")
    print(f"Total tracked: {len(covered)}")
    print(f"Total files: {len(all_files)}")
    print(f"Untracked files: {len(uncovered)}")
    
    if uncovered:
        print("\nUntracked files (may not trigger rebuilds):")
        by_ext = {}
        for f in uncovered:
            ext = f.suffix or "(no extension)"
            by_ext.setdefault(ext, []).append(f)
        
        for ext, files in sorted(by_ext.items()):
            print(f"\n  {ext}: {len(files)} files")
            for f in sorted(files)[:3]:
                print(f"    {f.relative_to(raw_dir)}")
            if len(files) > 3:
                print(f"    ... and {len(files)-3} more")
    
    print("\n" + "="*80)
    print("RECOMMENDATION")
    print("="*80)
    
    if len(uncovered) == 0:
        print("✓ All files are tracked by Snakemake rules")
    else:
        untracked_pct = 100 * len(uncovered) / len(all_files)
        print(f"⚠ {len(uncovered)} files ({untracked_pct:.1f}%) are not tracked")
        print("  These files won't trigger rebuilds when modified")
        print("  Check if any are important (e.g., images, PDFs with data)")

if __name__ == "__main__":
    main()