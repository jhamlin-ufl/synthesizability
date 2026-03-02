# src/synthesizability/oqmd.py
"""
OQMD database interface for querying hull distances and extracting structures.
"""
import subprocess
import re
from collections import defaultdict
from itertools import combinations
from typing import Optional, Dict, List
import numpy as np
from pymatgen.core import Structure, Lattice

OQMD_VERSION_REQUIRED = "1.7"
OQMD_DOWNLOAD_URL = "https://oqmd.org/download/"


def run_mysql_query(query: str) -> Optional[str]:
    """
    Execute a MySQL query on qmdb database.

    Args:
        query: SQL query string

    Returns:
        Query output as string, or None if failed
    """
    try:
        cmd = ['mysql', '-e', f'USE qmdb; {query}']
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        return None


def parse_formula_to_oqmd(formula: str) -> str:
    """
    Convert formula to OQMD format with alphabetically ordered elements.

    Args:
        formula: Chemical formula like 'MoNb2Ta1' or 'Al4FeCo3'

    Returns:
        OQMD formatted formula like 'Mo1 Nb2 Ta1' or 'Al4 Co3 Fe1'

    Examples:
        >>> parse_formula_to_oqmd('MoNb2Ta1')
        'Mo1 Nb2 Ta1'
        >>> parse_formula_to_oqmd('Al4FeCo3')
        'Al4 Co3 Fe1'
    """
    # Parse element-number pairs
    pattern = r'([A-Z][a-z]?)(\d*\.?\d*)'
    matches = re.findall(pattern, formula)

    # Build dict of element counts
    elem_dict = defaultdict(str)
    for element, count in matches:
        if element:
            count = count if count else '1'
            elem_dict[element] = count

    # Sort alphabetically and join
    parts = [f"{elem}{elem_dict[elem]}" for elem in sorted(elem_dict.keys())]
    return ' '.join(parts)


def parse_elements_from_formula(formula: str) -> list[str]:
    """
    Extract sorted list of unique elements from a formula string.

    Args:
        formula: Chemical formula like 'CoSi3Y' or 'MoNb2Ta'

    Returns:
        Sorted list of element symbols, e.g. ['Co', 'Si', 'Y']
    """
    pattern = r'([A-Z][a-z]?)(\d*\.?\d*)'
    matches = re.findall(pattern, formula)
    elements = sorted(set(elem for elem, _ in matches if elem))
    return elements


def query_formation_energies(oqmd_formula: str) -> List[Dict]:
    """
    Query OQMD database for formation energies of a composition.

    Args:
        oqmd_formula: Formula in OQMD format (e.g., 'Mo1 Nb2 Ta1')

    Returns:
        List of dicts with keys: composition_id, delta_e, stability, entry_id
    """
    query = f"""
    SELECT fe.composition_id, fe.delta_e, fe.stability, fe.entry_id
    FROM formation_energies fe
    WHERE fe.composition_id = '{oqmd_formula}'
    AND fe.fit_id = 'standard'
    """

    result = run_mysql_query(query)
    if not result:
        return []

    lines = result.strip().split('\n')
    if len(lines) <= 1:  # Only header or empty
        return []

    data = []
    for line in lines[1:]:  # Skip header
        parts = line.split('\t')
        if len(parts) == 4:
            data.append({
                'composition_id': parts[0],
                'delta_e': float(parts[1]) if parts[1] != 'NULL' else None,
                'stability': float(parts[2]) if parts[2] != 'NULL' else None,
                'entry_id': int(parts[3])
            })
    return data



def query_ternary_system(elements: list[str]) -> list[dict]:
    """
    Query all OQMD entries whose element set is a subset of the given elements.

    Covers all unaries, binaries, and ternaries (or fewer) within the chemical
    space defined by elements. Uses a three-step approach for performance:
    1. Fast subquery to get valid entry_ids via entries_element_set (~1.3s)
    2. Fast IN-clause lookup for formation energies (~0.02s)
    3. Fast IN-clause lookup for ICSD flags (~0.01s)

    Args:
        elements: Sorted list of element symbols, e.g. ['Co', 'Si', 'Y']

    Returns:
        List of dicts with keys: entry_id, composition_id, delta_e, stability,
        icsd (bool). One row per entry (all polymorphs included). Sorted by
        stability ascending (most stable first).
    """
    in_clause = ', '.join(f"'{el}'" for el in elements)

    # Step 1: get all entry_ids whose element set is a subset of our elements
    id_result = run_mysql_query(f"""
    SELECT entry_id FROM entries_element_set
    GROUP BY entry_id
    HAVING SUM(element_id NOT IN ({in_clause})) = 0
    """)
    if not id_result:
        return []

    entry_ids = [l.strip() for l in id_result.strip().split('\n')[1:] if l.strip()]
    if not entry_ids:
        return []

    ids_str = ', '.join(entry_ids)

    # Step 2: formation energies for those entry_ids
    fe_result = run_mysql_query(f"""
    SELECT e.id, e.composition_id, fe.delta_e, fe.stability
    FROM entries e
    JOIN formation_energies fe ON fe.entry_id = e.id
    WHERE fe.entry_id IN ({ids_str})
    AND fe.fit_id = 'standard'
    ORDER BY fe.stability ASC
    """)
    if not fe_result:
        return []

    fe_lines = fe_result.strip().split('\n')
    if len(fe_lines) <= 1:
        return []

    # Step 3: ICSD flags for those entry_ids
    icsd_result = run_mysql_query(f"""
    SELECT DISTINCT emd.entry_id
    FROM entries_meta_data emd
    JOIN meta_data md ON emd.metadata_id = md.id
    WHERE md.type='keyword' AND md.value='icsd'
    AND emd.entry_id IN ({ids_str})
    """)
    icsd_ids = set()
    if icsd_result:
        icsd_ids = set(
            int(l.strip()) for l in icsd_result.strip().split('\n')[1:]
            if l.strip().isdigit()
        )

    # Assemble results
    entries = []
    for line in fe_lines[1:]:
        parts = line.split('\t')
        if len(parts) == 4:
            eid = int(parts[0])
            entries.append({
                'entry_id': eid,
                'composition_id': parts[1],
                'delta_e': float(parts[2]) if parts[2] != 'NULL' else None,
                'stability': float(parts[3]) if parts[3] != 'NULL' else None,
                'icsd': eid in icsd_ids,
            })
    return entries

def query_exact_space(elements: list[str]) -> list[dict]:
    """
    Query all OQMD entries whose element set matches exactly the given elements.

    Unlike query_ternary_system, this returns only entries whose composition
    contains exactly these elements — no subsets. Use this to populate
    per-space CIF directories without duplication across subspaces.

    Uses the same three-step approach as query_ternary_system for performance.

    Args:
        elements: Sorted list of element symbols, e.g. ['Al', 'Co', 'Fe']

    Returns:
        List of dicts with keys: entry_id, composition_id, delta_e, stability,
        icsd (bool). One row per entry (all polymorphs included). Sorted by
        stability ascending (most stable first).
    """
    in_clause = ', '.join(f"'{el}'" for el in elements)
    n = len(elements)

    # Step 1: get entry_ids with exactly this element set
    id_result = run_mysql_query(f"""
    SELECT entry_id FROM entries_element_set
    GROUP BY entry_id
    HAVING SUM(element_id NOT IN ({in_clause})) = 0
    AND COUNT(DISTINCT element_id) = {n}
    """)
    if not id_result:
        return []

    entry_ids = [l.strip() for l in id_result.strip().split('\n')[1:] if l.strip()]
    if not entry_ids:
        return []

    ids_str = ', '.join(entry_ids)

    # Step 2: formation energies for those entry_ids
    fe_result = run_mysql_query(f"""
    SELECT e.id, e.composition_id, fe.delta_e, fe.stability
    FROM entries e
    JOIN formation_energies fe ON fe.entry_id = e.id
    WHERE fe.entry_id IN ({ids_str})
    AND fe.fit_id = 'standard'
    ORDER BY fe.stability ASC
    """)
    if not fe_result:
        return []

    fe_lines = fe_result.strip().split('\n')
    if len(fe_lines) <= 1:
        return []

    # Step 3: ICSD flags for those entry_ids
    icsd_result = run_mysql_query(f"""
    SELECT DISTINCT emd.entry_id
    FROM entries_meta_data emd
    JOIN meta_data md ON emd.metadata_id = md.id
    WHERE md.type='keyword' AND md.value='icsd'
    AND emd.entry_id IN ({ids_str})
    """)
    icsd_ids = set()
    if icsd_result:
        icsd_ids = set(
            int(l.strip()) for l in icsd_result.strip().split('\n')[1:]
            if l.strip().isdigit()
        )

    entries = []
    for line in fe_lines[1:]:
        parts = line.split('\t')
        if len(parts) == 4:
            eid = int(parts[0])
            entries.append({
                'entry_id': eid,
                'composition_id': parts[1],
                'delta_e': float(parts[2]) if parts[2] != 'NULL' else None,
                'stability': float(parts[3]) if parts[3] != 'NULL' else None,
                'icsd': eid in icsd_ids,
            })
    return entries

def make_cif_filename(composition_id: str, entry_id: int, stability: Optional[float]) -> str:
    """
    Build a CIF filename encoding composition, entry ID, and stability.

    Stability is stored in meV/atom, rounded to integer, to avoid extra
    decimal points in the filename. NULL stability is indicated explicitly.

    Args:
        composition_id: OQMD composition string, e.g. 'Co2 Si2 Y1'
        entry_id: OQMD entry ID
        stability: Stability in eV/atom, or None

    Returns:
        Filename string, e.g. 'Co2Si2Y1_647098_stab-92meV.cif'
    """
    # Compact composition: remove spaces
    comp_compact = composition_id.replace(' ', '')

    if stability is not None:
        stab_meV = round(stability * 1000)
        stab_str = f"stab{stab_meV:+d}meV"
    else:
        stab_str = "stabUNKNOWN"

    return f"{comp_compact}_{entry_id}_{stab_str}.cif"


def get_structure_id_for_entry(entry_id: int) -> Optional[int]:
    """
    Get the output structure ID for a given entry_id.

    Args:
        entry_id: OQMD entry ID

    Returns:
        Structure ID (output_id from calculations), or None if not found
    """
    query = f"""
    SELECT c.output_id
    FROM formation_energies fe
    JOIN calculations c ON fe.calculation_id = c.id
    WHERE fe.entry_id = {entry_id} AND fe.fit_id = 'standard'
    """

    result = run_mysql_query(query)
    if not result:
        return None

    lines = result.strip().split('\n')
    if len(lines) > 1:
        return int(lines[1])
    return None


def get_structure_from_db(entry_id: int) -> Optional[Structure]:
    """
    Extract pymatgen Structure object from OQMD database.

    Args:
        entry_id: OQMD entry ID

    Returns:
        pymatgen Structure object, or None if extraction fails
    """
    # Get structure_id
    structure_id = get_structure_id_for_entry(entry_id)
    if structure_id is None:
        return None

    # Get lattice vectors
    lattice_query = f"""
    SELECT x1, x2, x3, y1, y2, y3, z1, z2, z3
    FROM structures
    WHERE id = {structure_id}
    """

    lattice_result = run_mysql_query(lattice_query)
    if not lattice_result:
        return None

    lines = lattice_result.strip().split('\n')
    if len(lines) <= 1:
        return None

    lat_parts = lines[1].split('\t')
    if len(lat_parts) != 9:
        return None

    # Build lattice matrix [[x1,x2,x3], [y1,y2,y3], [z1,z2,z3]]
    lattice_matrix = [
        [float(lat_parts[0]), float(lat_parts[1]), float(lat_parts[2])],
        [float(lat_parts[3]), float(lat_parts[4]), float(lat_parts[5])],
        [float(lat_parts[6]), float(lat_parts[7]), float(lat_parts[8])]
    ]
    lattice = Lattice(lattice_matrix)

    # Get atomic positions
    atoms_query = f"""
    SELECT element_id, x, y, z
    FROM atoms
    WHERE structure_id = {structure_id}
    ORDER BY id
    """

    atoms_result = run_mysql_query(atoms_query)
    if not atoms_result:
        return None

    lines = atoms_result.strip().split('\n')
    if len(lines) <= 1:
        return None

    species = []
    coords = []
    for line in lines[1:]:
        parts = line.split('\t')
        if len(parts) == 4:
            species.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    if not species:
        return None

    structure = Structure(lattice, species, coords)
    return structure


def check_database_exists() -> bool:
    """Check if qmdb database exists."""
    try:
        cmd = ['mysql', '-e', 'SHOW DATABASES LIKE "qmdb";']
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return 'qmdb' in result.stdout
    except:
        return False


def get_database_entry_count() -> Optional[int]:
    """Get number of entries in OQMD database."""
    result = run_mysql_query("SELECT COUNT(*) FROM entries;")
    if result:
        lines = result.split('\n')
        if len(lines) > 1:
            return int(lines[1])
    return None
