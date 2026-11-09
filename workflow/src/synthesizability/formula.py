"""Formula presentation: element order by Allen electronegativity, subscripts."""
import re

# Allen electronegativity (configuration energy).  Under this scale Si, Ge and Sn
# sit above the late transition metals, giving the conventional order for these
# intermetallics (ZrPt5Si, not ZrSiPt5).  La and Gd are estimates (Allen omits
# the f block); both are the most electropositive element wherever they occur.
ALLEN = {
    'Al': 1.613, 'Au': 1.92, 'B': 2.051, 'C': 2.544, 'Co': 1.84, 'Cr': 1.65,
    'Cu': 1.85, 'Fe': 1.80, 'Gd': 1.20, 'Ge': 1.994, 'Hf': 1.16, 'Ir': 1.68,
    'La': 1.10, 'Mn': 1.75, 'Mo': 1.47, 'Nb': 1.41, 'Ni': 1.88, 'Pd': 1.58,
    'Pt': 1.72, 'Ru': 1.54, 'Sc': 1.19, 'Si': 1.916, 'Sn': 1.824, 'Ti': 1.38,
    'Y': 1.12, 'Zr': 1.32,
}


def parse(formula):
    """'ZrSiPt5' -> [('Zr', ''), ('Si', ''), ('Pt', '5')]."""
    return re.findall(r'([A-Z][a-z]?)(\d*)', formula)


def ordered(formula):
    """'ZrSiPt5' -> 'ZrPt5Si' (increasing electronegativity)."""
    return ''.join(el + n for el, n in sorted(parse(formula), key=lambda p: (ALLEN[p[0]], p[0])))


def mathtext(formula):
    """'ZrSiPt5' -> '$\\mathrm{ZrPt_{5}Si}$' for matplotlib labels."""
    return '$\\mathrm{' + ''.join(f'{el}_{{{n}}}' if n else el for el, n in parse(ordered(formula))) + '}$'

def parse_formula(formula):
    """'ZrSiPt5' -> [('Zr', ''), ('Si', ''), ('Pt', '5')]."""
    return [(el, n) for el, n in re.findall(r'([A-Z][a-z]?)(\d*)', formula) if el]


def sort_by_electronegativity(formula):
    """'ZrSiPt5' -> 'ZrPt5Si': elements in order of increasing electronegativity."""
    parts = parse_formula(formula)
    missing = [el for el, _ in parts if el not in ALLEN_ELECTRONEGATIVITY]
    if missing:
        raise KeyError(f'no Allen electronegativity for {missing} '
                       f'(formula {formula}); add it to ALLEN_ELECTRONEGATIVITY')
    parts.sort(key=lambda p: (ALLEN_ELECTRONEGATIVITY[p[0]], p[0]))
    return ''.join(el + n for el, n in parts)


def formula_to_mathtext(formula):
    """'ZrSiPt5' -> '$\\mathrm{ZrPt_{5}Si}$': reordered, with subscripts."""
    return ('$\\mathrm{'
            + ''.join(f'{el}_{{{n}}}' if n else el
                      for el, n in parse_formula(sort_by_electronegativity(formula)))
            + '}$')
