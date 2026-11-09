"""Mixture hypotheses from convex-hull geometry.

A sample of nominal composition c that did not form the target decomposes,
at equilibrium, into the phases whose simplex on the hull contains c.  The
energy of any mixture of phases is the barycentric average of their energies,
so every three-phase triangle (or two-phase segment, or single phase) that
contains c is a hypothesis with a well-defined energy excess over the hull,
and the barycentric coordinates are the lever-rule phase fractions per atom.
The excess is used as a prior for ranking, not as a hard filter.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import linprog

TOL = 1e-9
MIN_FRACTION = 1e-3   # a vertex with less than this is not a distinct phase
# Real buttons are a little off the nominal composition (weighing, mass loss),
# so simplices whose barycentric coordinates are slightly negative at c are
# kept too: they are the "target plus a small amount of a neighbour" cases
# that an exact-composition test can never propose.  Negative coordinates are
# clipped to IMPURITY_FLOOR and the fractions renormalised.
OFF_NOMINAL = 0.15
IMPURITY_FLOOR = 0.03


def hull_at(X, E, c):
    """Hull energy per atom at c and the LP weights (vertex solution).

    Only two composition coordinates are constrained (the third follows from
    the sum rule), so rounded fractions cannot make the system inconsistent.
    """
    n = len(E)
    X = np.asarray(X, float)
    c = np.asarray(c, float)
    res = linprog(E, A_eq=np.vstack([X[:, :2].T, np.ones(n)]), b_eq=np.append(c[:2], 1.0),
                  bounds=(0, None), method='highs')
    if not res.success:
        raise RuntimeError(f'hull LP failed: {res.message}')
    return float(res.fun), res.x


def composition_points(phases):
    """Unique composition points over all sources.

    Returns (X (m, 3), members: list of index lists into ``phases``).
    """
    groups = {}
    for i, p in enumerate(phases):
        groups.setdefault(tuple(np.round(p.x, 6)), []).append(i)
    X = np.array(list(groups.keys()))
    return X, list(groups.values())


def source_energies(phases, members, source):
    """Lowest formation energy of ``source`` at each composition point (nan
    where the source has no entry).  Sources sharing an energy scale
    (``library.SOURCE_GROUPS``) are pooled."""
    from .library import SOURCE_GROUPS
    group = SOURCE_GROUPS.get(source, (source,))
    E = np.full(len(members), np.nan)
    for n, idx in enumerate(members):
        es = [phases[i].e_form for i in idx
              if phases[i].source in group and phases[i].e_form is not None]
        if es:
            E[n] = min(es)
    return E


def enumerate_all(phases, c, sources=('oqmd', 'alex', 'mp'), max_excess=None):
    """Mixture hypotheses from every source's hull, merged on composition.

    Each source's formation energies are only comparable among themselves, so
    the hull geometry is evaluated per source over the composition points that
    source covers, then the hypotheses are merged; a mixture found by several
    sources keeps its smallest excess and the source that gave it.  Point
    indices refer to ``composition_points(phases)``.
    """
    X, members = composition_points(phases)
    frames = []
    hulls = {}
    for source in sources:
        E = source_energies(phases, members, source)
        have = np.where(np.isfinite(E))[0]
        if have.size < 3:
            continue
        try:
            df = enumerate_mixtures(X[have], E[have], c, max_excess)
        except RuntimeError:
            continue
        hulls[source] = df.attrs['hull']
        for col in ('i', 'j', 'k'):
            df[col] = np.where(df[col] >= 0, have[np.clip(df[col].to_numpy(), 0, None)], -1)
        df['prior_source'] = source
        frames.append(df)
    if not frames:
        raise RuntimeError('no source has enough energies for a hull')
    allm = pd.concat(frames, ignore_index=True)
    # Canonical vertex order so the same triple from two sources merges
    pts = allm[['i', 'j', 'k']].to_numpy()
    fr = allm[['f1', 'f2', 'f3']].to_numpy()
    order = np.argsort(np.where(pts < 0, 10 ** 9, pts), axis=1)
    pts = np.take_along_axis(pts, order, axis=1)
    fr = np.take_along_axis(fr, order, axis=1)
    allm[['i', 'j', 'k']] = pts
    allm[['f1', 'f2', 'f3']] = fr
    merged = (allm.sort_values('excess').groupby(['i', 'j', 'k'], as_index=False).first()
              .sort_values('excess').reset_index(drop=True))
    merged.attrs['hulls'] = hulls
    return X, members, merged


def enumerate_mixtures(X, E, c, max_excess=None):
    """All simplices of composition points that contain c.

    Returns a DataFrame with the point indices (``i``, ``j``, ``k``; -1 for
    absent vertices), lever-rule fractions ``f1 f2 f3``, mixture energy
    ``e_mix`` and ``excess`` over the hull at c (eV/atom).
    """
    hull, _ = hull_at(X, E, c)
    m = len(E)
    rows = []

    # Single phases at c
    at_c = np.where(np.all(np.abs(X - c) < 1e-6, axis=1))[0]
    for i in at_c:
        rows.append((i, -1, -1, 1.0, 0.0, 0.0, E[i], 0.0))

    # Two-phase segments through c (c strictly between the end points)
    for i in range(m):
        d = X[i] - c
        if np.all(np.abs(d) < 1e-6):
            continue
        for j in range(i + 1, m):
            dj = X[j] - c
            if np.all(np.abs(dj) < 1e-6):
                continue
            # collinear and on opposite sides: dj = -t d with t > 0
            cross = d[0] * dj[1] - d[1] * dj[0]
            if abs(cross) > 1e-9:
                continue
            k = np.argmax(np.abs(d))
            t = -dj[k] / d[k]
            f1 = t / (1 + t)
            off = -min(f1, 1 - f1)
            if off > OFF_NOMINAL:
                continue
            if off > -MIN_FRACTION:
                # c sits at (or just outside) an end point: impurity case
                f1 = max(f1, IMPURITY_FLOOR) if f1 < 0.5 else min(f1, 1 - IMPURITY_FLOOR)
            rows.append((i, j, -1, f1, 1 - f1, 0.0, f1 * E[i] + (1 - f1) * E[j], max(off, 0.0)))

    # Triangles: for each (i, j) solve for every k at once
    idx = np.arange(m)
    for i in range(m):
        for j in range(i + 1, m):
            k = idx[j + 1:]
            if k.size == 0:
                continue
            xk = X[k][:, :2]
            M = np.stack([X[i][:2] - xk, X[j][:2] - xk], axis=2)
            det = np.linalg.det(M)
            ok = np.abs(det) > 1e-12
            if not ok.any():
                continue
            l12 = np.full((k.size, 2), np.nan)
            l12[ok] = np.linalg.solve(M[ok], (c[:2] - xk)[ok][..., None])[..., 0]
            l3 = 1 - l12.sum(1)
            L = np.column_stack([l12, l3])
            lmin = np.nanmin(L, axis=1)
            inside = ok & (lmin > -OFF_NOMINAL)
            if not inside.any():
                continue
            kk = k[inside]
            Lc = np.clip(L[inside], IMPURITY_FLOOR, None)
            Lc = Lc / Lc.sum(1, keepdims=True)
            # Simplices where c is on an edge or vertex duplicate lower-order
            # hypotheses unless they add an impurity; keep them flagged
            off = np.clip(-lmin[inside], 0.0, None)
            f1, f2, f3 = Lc[:, 0], Lc[:, 1], Lc[:, 2]
            e_mix = f1 * E[i] + f2 * E[j] + f3 * E[kk]
            for row in zip(np.full(kk.size, i), np.full(kk.size, j), kk, f1, f2, f3, e_mix, off):
                rows.append(row)

    df = pd.DataFrame(rows, columns=['i', 'j', 'k', 'f1', 'f2', 'f3', 'e_mix', 'off_nominal'])
    # Impurity hypotheses sit slightly off c, so their energy can fall below
    # the hull at c; the prior for a known-stable target plus impurities is 0
    df['excess'] = np.maximum(df['e_mix'] - hull, 0.0)
    if max_excess is not None:
        df = df[df['excess'] <= max_excess + TOL]
    df = df.sort_values('excess').reset_index(drop=True)
    df.attrs['hull'] = hull
    return df


def hull_excluding_target(X, E, c):
    """Hull energy at c with the points at c removed (decomposition energy)."""
    keep = ~np.all(np.abs(X - c) < 1e-6, axis=1)
    if keep.sum() < 3:
        return None
    hull, _ = hull_at(X[keep], E[keep], c)
    return hull


def atom_to_weight_fractions(fracs, comps):
    """Lever-rule fractions per atom -> weight fractions."""
    masses = np.array([c.weight / c.num_atoms for c in comps])
    w = np.asarray(fracs) * masses
    return w / w.sum()
