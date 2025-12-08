#!/usr/bin/env python3
"""
SuperConductor Element Finder

Searches the SuperCon database for superconducting materials containing
only specified elements.

Data source: MDR SuperCon Datasheet Ver.240322
    Publisher: National Institute for Materials Science (NIMS)
    DOI: https://doi.org/10.48505/nims.3837
    License: CC BY 4.0

Citation:
    Materials Database Group. MDR SuperCon Datasheet Ver.240322.
    https://doi.org/10.48505/nims.3837

Usage:
    python find_relevant_SCs.py Mo Nb Ta
    python find_relevant_SCs.py --include-oxygen La Ba Cu
"""

import pandas as pd
import argparse
import re
from pathlib import Path
from pymatgen.core import Composition
import urllib.request
import urllib.parse
import zipfile
import sys

__version__ = "1.0.0"
__author__ = "James J. Hamlin"

# Paths relative to repository root
REPO_ROOT = Path(__file__).parent.parent
DATA_DIR = REPO_ROOT / "data/external/supercon"
DATA_FILE = DATA_DIR / "primary.tsv"
DATA_URL = "https://mdr.nims.go.jp/datasets/5d8000f3-a8cd-4ad5-bcdb-2447a5166839.zip"


def ensure_data_available():
    """
    Download SuperCon database if not present.
    
    Downloads and extracts the MDR SuperCon Datasheet from NIMS.
    Data is stored in data/external/supercon/ relative to repository root.
    """
    if DATA_FILE.exists():
        return
    
    print("SuperCon database not found.", file=sys.stderr)
    print(f"Downloading from {DATA_URL}...", file=sys.stderr)
    
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    
    zip_path = DATA_DIR / "supercon_data.zip"
    
    try:
        urllib.request.urlretrieve(DATA_URL, zip_path)
    except Exception as e:
        print(f"Error downloading data: {e}", file=sys.stderr)
        sys.exit(1)
    
    print("Extracting files...", file=sys.stderr)
    
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(DATA_DIR)
        zip_path.unlink()
    except Exception as e:
        print(f"Error extracting data: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Database downloaded to {DATA_DIR}\n", file=sys.stderr)


def extract_doi(journal_str):
    """
    Extract or construct DOI from journal citation.
    
    Attempts to construct DOIs for common physics journals where the
    format is predictable (Physical Review, J.Phys.Soc.Jpn).
    
    Args:
        journal_str: Journal citation string from SuperCon database
        
    Returns:
        DOI string if constructable, None otherwise
    """
    if pd.isna(journal_str):
        return None
    
    cite = str(journal_str).strip()
    
    # Check if DOI already present
    doi_match = re.search(r'10\.\d{4,}/[^\s,]+', cite)
    if doi_match:
        return doi_match.group(0)
    
    # Physical Review B: "Phys.Rev.B,57(1998)7491"
    prb = re.search(r'Phys\.?Rev\.?B[,\s]+(\d+)\s*\((\d{4})\)\s*(\d+)', cite, re.I)
    if prb:
        vol, year, page = prb.groups()
        return f"10.1103/PhysRevB.{vol}.{page}"
    
    # Physical Review Letters: "Phys.Rev.Lett.,77(1998)163"
    prl = re.search(r'Phys\.?Rev\.?Lett\.?[,\s]+(\d+)\s*\((\d{4})\)\s*(\d+)', cite, re.I)
    if prl:
        vol, year, page = prl.groups()
        return f"10.1103/PhysRevLett.{vol}.{page}"
    
    # Physical Review (pre-1970): "PHYS.REV.,165(1968)533"
    pr = re.search(r'PHYS\.?REV\.[,\s]+(\d+)\s*\((\d{4})\)\s*([A-Z]?\d+)', cite, re.I)
    if pr:
        vol, year, page = pr.groups()
        if int(year) >= 1970:
            return f"10.1103/PhysRevB.{vol}.{page}"
        else:
            return f"10.1103/PhysRev.{vol}.{page}"
    
    # J.Phys.Soc.Japan: "J.PHYS.SOC.JAPAN,56(1987)3805"
    jpsj = re.search(r'J\.PHYS\.SOC\.JAPAN[,\s]+(\d+)\s*\((\d{4})\)\s*(\d+)', cite, re.I)
    if jpsj:
        vol, year, page = jpsj.groups()
        return f"10.1143/JPSJ.{vol}.{page}"
    
    return None


def get_reference_link(journal_str, doi):
    """
    Generate reference link for citation.
    
    Returns DOI link if available, otherwise generates a Google Scholar
    search link to help locate the paper.
    
    Args:
        journal_str: Journal citation string
        doi: DOI string or None
        
    Returns:
        URL string for reference lookup
    """
    if doi:
        return f"https://doi.org/{doi}"
    
    if pd.isna(journal_str):
        return None
    
    # Use Google Scholar search for citations without DOI
    query = urllib.parse.quote(str(journal_str))
    return f"https://scholar.google.com/scholar?q={query}"


def extract_elements_from_formula(formula):
    """
    Extract element symbols from chemical formula.
    
    Uses pymatgen to parse chemical formulas and extract elements,
    excluding oxygen by default.
    
    Args:
        formula: Chemical formula string (e.g., "Ba0.2La1.8Cu1O4-Y")
        
    Returns:
        Set of element symbols (strings)
    """
    try:
        # Remove trailing markers like -Y, -Z
        formula_clean = re.sub(r'-[A-Z]$', '', formula)
        comp = Composition(formula_clean)
        elements = {elem.symbol for elem in comp.elements}
        elements.discard('O')
        return elements
    except Exception:
        return set()


def find_matching_superconductors(target_elements, include_oxygen=False):
    """
    Search SuperCon database for materials matching element criteria.
    
    Args:
        target_elements: List of element symbols to search for
        include_oxygen: If True, include oxygen-containing compounds
        
    Returns:
        List of dictionaries containing matched superconductor data
    """
    ensure_data_available()
    
    # Load SuperCon primary data
    df = pd.read_csv(DATA_FILE, sep='\t', skiprows=2)
    df.columns = ['num', 'name', 'element', 'str3', 'utc', 'tc', 'journal']
    
    target_set = set(target_elements)
    matches = []
    
    for idx, row in df.iterrows():
        compound_elements = extract_elements_from_formula(row['element'])
        
        # Skip oxygen-containing compounds unless requested
        if not include_oxygen and 'O' in compound_elements:
            continue
        
        # Check if compound elements are subset of target elements
        if compound_elements and compound_elements.issubset(target_set):
            matches.append({
                'compound': row['name'],
                'formula': row['element'],
                'elements': sorted(compound_elements),
                'tc': row['tc'],
                'journal': row['journal'],
                'doi': extract_doi(row['journal'])
            })
    
    return matches


def format_output(matches):
    """
    Format search results for display.
    
    Args:
        matches: List of match dictionaries from find_matching_superconductors
    """
    print(f"\nFound {len(matches)} superconducting materials")
    print("Note: Tc values may include measurements under pressure.")
    print("      Consult original references for experimental conditions.\n")
    print("-" * 140)
    print(f"{'Material':<25} {'Formula':<30} {'Elements':<15} {'Tc (K)':<8} {'Reference'}")
    print("-" * 140)
    
    # Sort by Tc (highest first)
    matches.sort(key=lambda x: x['tc'] if pd.notna(x['tc']) else -999, reverse=True)
    
    for match in matches:
        elements_str = '+'.join(match['elements'])
        tc_str = f"{match['tc']:.1f}" if pd.notna(match['tc']) else "N/A"
        
        # Use formula as name if common name is missing
        material_name = match['compound'] if pd.notna(match['compound']) else match['formula']
        
        journal_str = str(match['journal']) if pd.notna(match['journal']) else "N/A"
        
        # Add reference link (DOI or Google Scholar search)
        link = get_reference_link(match['journal'], match['doi'])
        if link:
            journal_str = f"{journal_str} {link}"
        
        print(f"{material_name:<25} {match['formula']:<30} {elements_str:<15} "
              f"{tc_str:<8} {journal_str}")
    
    print(f"\nTotal: {len(matches)} materials\n")


def main():
    """Main entry point for command-line interface."""
    parser = argparse.ArgumentParser(
        description='Search SuperCon database for superconductors by element composition',
        epilog='Data source: MDR SuperCon Datasheet (NIMS) https://doi.org/10.48505/nims.3837',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        'elements',
        nargs='+',
        help='Element symbols to search for (e.g., Mo Nb Ta)'
    )
    
    parser.add_argument(
        '--include-oxygen',
        action='store_true',
        help='Include oxygen-containing compounds in results'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {__version__}'
    )
    
    args = parser.parse_args()
    
    # Validate element symbols
    for elem in args.elements:
        if not elem[0].isupper() or (len(elem) > 1 and not elem[1:].islower()):
            print(f"Warning: '{elem}' may not be a valid element symbol", file=sys.stderr)
    
    print(f"Searching for superconductors containing: {', '.join(sorted(args.elements))}")
    if not args.include_oxygen:
        print("(excluding oxygen-containing compounds)")
    
    matches = find_matching_superconductors(args.elements, args.include_oxygen)
    format_output(matches)


if __name__ == "__main__":
    main()