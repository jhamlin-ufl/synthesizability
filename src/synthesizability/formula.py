"""
Chemical formula analysis and property calculation.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from pymatgen.core import Composition, Element


# Load reference data at module scope
_DATA_DIR = Path(__file__).parent.parent.parent / "data" / "external" / "reference"

_price_table = pd.read_csv(_DATA_DIR / "element_prices.csv").set_index("Symb")
_vapor_pressure_data = pd.read_csv(_DATA_DIR / "element_vapor_pressures.csv").set_index("symbol")


def extract_formula_from_sample_id(sample_id: str) -> str:
    """
    Extract chemical formula from sample ID.
    
    Example: '0447_HM_MoTiTa2' -> 'MoTiTa2'
    """
    parts = sample_id.split('_')
    if len(parts) >= 3:
        formula = '_'.join(parts[2:])
        if '(' in formula:
            formula = formula.split('(')[0]
        return formula.strip()
    return ""


def calculate_price_per_gram(composition: Composition) -> float:
    """Calculate estimated price per gram of composition."""
    try:
        elements = [elem.symbol for elem in composition.elements]
        element_prices = _price_table.loc[elements]["$/gram"]
        fractions = [composition.get_atomic_fraction(elem) for elem in composition.elements]
        return np.multiply(element_prices, fractions).sum()
    except KeyError:
        return np.nan


def is_arc_meltable(composition: Composition) -> bool:
    """
    Determine if composition can be arc melted.
    
    Criteria:
    - All elements solid at room temp (MP > 300K)
    - All elements have vapor pressure < 1 torr at highest melting point
    """
    melting_points = np.array([elem.melting_point for elem in composition.elements])
    
    # Check if solid at room temp
    if np.any(melting_points < 300):
        return False
    
    highest_mp = np.max(melting_points)
    
    try:
        elements = [elem.symbol for elem in composition.elements]
        T_for_1torr = np.array([_vapor_pressure_data.loc[elem]["1torr"] for elem in elements])
        return np.all(T_for_1torr > highest_mp)
    except KeyError:
        return False


def enrich_with_formula_properties(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add formula-derived properties to dataframe.
    
    Adds columns: formula, price_per_gram, arc_meltable, disorder_probability
    """
    formulas = []
    prices = []
    arc_meltable_flags = []
    
    for sample_id in df['sample_id']:
        formula = extract_formula_from_sample_id(sample_id)
        formulas.append(formula)
        
        if formula:
            try:
                comp = Composition(formula)
                prices.append(calculate_price_per_gram(comp))
                arc_meltable_flags.append(is_arc_meltable(comp))
            except:
                prices.append(np.nan)
                arc_meltable_flags.append(False)
        else:
            prices.append(np.nan)
            arc_meltable_flags.append(False)
    
    df['formula'] = formulas
    df['price_per_gram'] = prices
    df['arc_meltable'] = arc_meltable_flags
    
    # Add disorder probabilities from cache
    df = add_disorder_probabilities(df)
    
    return df

def add_disorder_probabilities(df: pd.DataFrame, cache_path: Path = None) -> pd.DataFrame:
    """
    Add disorder probabilities to dataframe from cache.
    
    Args:
        df: Dataframe with 'formula' column
        cache_path: Path to disorder cache CSV (default: data/processed/disorder_cache.csv)
        
    Returns:
        Dataframe with added 'disorder_probability' column
    """
    if cache_path is None:
        cache_path = Path(__file__).parent.parent.parent / 'data' / 'processed' / 'disorder_cache.csv'
    
    if not cache_path.exists():
        print(f"WARNING: Disorder cache not found at {cache_path}")
        print("         Run 'snakemake data/processed/disorder_cache.csv' to compute disorder")
        df['disorder_probability'] = np.nan
        return df
    
    # Load cache
    cache_df = pd.read_csv(cache_path)
    
    # Merge with dataframe
    df = df.merge(cache_df, on='formula', how='left')
    
    # Report coverage
    n_with_disorder = df['disorder_probability'].notna().sum()
    print(f"Disorder probabilities: {n_with_disorder}/{len(df)} samples")
    
    return df
