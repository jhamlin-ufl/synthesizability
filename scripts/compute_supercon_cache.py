#!/usr/bin/env python3
"""
Compute SuperCon cache for all sample formulas.

For each sample formula in formulas.txt, queries the local SuperCon database
for all superconductors whose elements are a subset of the sample's elements.
Results are written as JSON files to results/supercon/<formula>.json.

Data source: MDR SuperCon Datasheet Ver.240322
    Publisher: National Institute for Materials Science (NIMS)
    DOI: https://doi.org/10.48505/nims.3837
"""

import json
import re
import sys
import urllib.parse
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition

REPO_ROOT = Path(__file__).parent.parent
DATA_FILE = REPO_ROOT / "data/external/supercon/primary.tsv"
FORMULAS_FILE = REPO_ROOT / "data/processed/formulas.txt"
OUTPUT_DIR = REPO_ROOT / "results/supercon"


def extract_elements(formula: str) -> set[str]:
    """Extract element symbols from a chemical formula, excluding oxygen."""
    try:
        formula_clean = re.sub(r'-[A-Z]$', '', formula)
        comp = Composition(formula_clean)
        elements = {elem.symbol for elem in comp.elements}
        elements.discard('O')
        return elements
    except Exception:
        return set()


def extract_doi(journal_str: str) -> str | None:
    """
    Extract or construct DOI from journal citation string.
    Handles Physical Review B/Letters, Physical Review, J.Phys.Soc.Japan.
    """
    if pd.isna(journal_str):
        return None

    cite = str(journal_str).strip()

    doi_match = re.search(r'10\.\d{4,}/[^\s,]+', cite)
    if doi_match:
        return doi_match.group(0)

    prb = re.search(r'Phys\.?Rev\.?B[,\s]+(\d+)\s*\((\d{4})\)\s*(\d+)', cite, re.I)
    if prb:
        vol, year, page = prb.groups()
        return f"10.1103/PhysRevB.{vol}.{page}"

    prl = re.search(r'Phys\.?Rev\.?Lett\.?[,\s]+(\d+)\s*\((\d{4})\)\s*(\d+)', cite, re.I)
    if prl:
        vol, year, page = prl.groups()
        return f"10.1103/PhysRevLett.{vol}.{page}"

    pr = re.search(r'PHYS\.?REV\.[,\s]+(\d+)\s*\((\d{4})\)\s*([A-Z]?\d+)', cite, re.I)
    if pr:
        vol, year, page = pr.groups()
        if int(year) >= 1970:
            return f"10.1103/PhysRevB.{vol}.{page}"
        else:
            return f"10.1103/PhysRev.{vol}.{page}"

    jpsj = re.search(r'J\.PHYS\.SOC\.JAPAN[,\s]+(\d+)\s*\((\d{4})\)\s*(\d+)', cite, re.I)
    if jpsj:
        vol, year, page = jpsj.groups()
        return f"10.1143/JPSJ.{vol}.{page}"

    return None


def get_reference_url(journal_str: str, doi: str | None) -> str | None:
    """Return DOI link if available, otherwise a Google Scholar search URL."""
    if doi:
        return f"https://doi.org/{doi}"
    if pd.isna(journal_str):
        return None
    query = urllib.parse.quote(str(journal_str))
    return f"https://scholar.google.com/scholar?q={query}"


def load_supercon() -> pd.DataFrame:
    """Load and preprocess the SuperCon TSV database."""
    df = pd.read_csv(DATA_FILE, sep='\t', skiprows=2)
    df.columns = ['num', 'name', 'element', 'str3', 'utc', 'tc', 'journal']
    df['_elements'] = df['element'].apply(extract_elements)
    # Drop rows where element parsing failed entirely
    df = df[df['_elements'].map(len) > 0].copy()
    print(f"Loaded SuperCon: {len(df)} entries after element parsing")
    return df


def query_for_elements(supercon_df: pd.DataFrame, target_elements: set[str]) -> list[dict]:
    """
    Return all SuperCon entries whose elements are a subset of target_elements,
    excluding oxygen-containing compounds. Sorted by Tc descending.
    """
    mask = supercon_df['_elements'].apply(
        lambda elems: bool(elems) and elems.issubset(target_elements)
    )
    hits = supercon_df[mask].copy()

    results = []
    for _, row in hits.iterrows():
        doi = extract_doi(row['journal'])
        results.append({
            'compound': row['name'] if pd.notna(row['name']) else None,
            'formula': row['element'],
            'tc': float(row['tc']) if pd.notna(row['tc']) else None,
            'doi': doi,
            'journal': str(row['journal']) if pd.notna(row['journal']) else None,
            'url': get_reference_url(row['journal'], doi),
        })

    results.sort(key=lambda x: x['tc'] if x['tc'] is not None else -999, reverse=True)
    return results


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    formulas = FORMULAS_FILE.read_text().strip().splitlines()
    print(f"Processing {len(formulas)} formulas")

    supercon_df = load_supercon()

    for i, formula in enumerate(formulas, 1):
        target_elements = extract_elements(formula)
        if not target_elements:
            print(f"  [{i}/{len(formulas)}] {formula}: could not parse elements, skipping")
            continue

        hits = query_for_elements(supercon_df, target_elements)
        out_path = OUTPUT_DIR / f"{formula}.json"
        out_path.write_text(json.dumps(hits, indent=2))
        print(f"  [{i}/{len(formulas)}] {formula}: {len(hits)} hits → {out_path.name}")

    print(f"\nDone. JSON files written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()