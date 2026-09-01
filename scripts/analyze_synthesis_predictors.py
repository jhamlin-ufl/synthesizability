#!/usr/bin/env python3
"""
How well do hull distance, disorder parameter, and dynamic stability predict
whether an arc-melted target forms in its predicted structure?

Operates on the 33 classified targets of Fig. disorder_vs_stability (one row per
target, so remakes are not double-counted).

Two outcomes are analysed:
  success  — the predicted structure was found (P)
  disorder — the sample formed a solid solution (SS)

For each predictor we report the discrimination (AUC via Mann-Whitney U, with a
bootstrap CI), a significance test, and leave-one-out cross-validated AUC for
single-feature and combined logistic models.  With n=33 and only 6 successes
this is an underpowered dataset; the cross-validated numbers are the honest
measure of predictive ability and the multivariate fit is indicative only.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).parent))
from plot_disorder_vs_stability import load_classified_samples, load_dynamic_stability

RNG = np.random.default_rng(0)
N_BOOT = 10000

PREDICTORS = {
    'disorder_probability': 'Disorder parameter',
    'oqmd_stability':       'Hull distance (eV/atom)',
    'dyn_unstable':         'Ground state dynamically unstable',
}


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def auc_and_test(x, y):
    """AUC for x separating y, plus Mann-Whitney U p-value.

    AUC > 0.5 means larger x goes with y == True.
    """
    a, b = np.asarray(x)[y], np.asarray(x)[~y]
    u = stats.mannwhitneyu(a, b, alternative='two-sided')
    return u.statistic / (len(a) * len(b)), u.pvalue


def auc_ci(x, y, n_boot=N_BOOT, alpha=0.05):
    """Percentile bootstrap CI for the AUC, resampling within each class."""
    x = np.asarray(x, dtype=float)
    idx_pos, idx_neg = np.flatnonzero(y), np.flatnonzero(~y)
    boots = []
    for _ in range(n_boot):
        p = RNG.choice(idx_pos, len(idx_pos), replace=True)
        n = RNG.choice(idx_neg, len(idx_neg), replace=True)
        a, b = x[p], x[n]
        # AUC as P(a > b) + 0.5 P(a == b)
        diff = a[:, None] - b[None, :]
        boots.append((diff > 0).mean() + 0.5 * (diff == 0).mean())
    return np.quantile(boots, [alpha / 2, 1 - alpha / 2])


def loo_auc(X, y):
    """Leave-one-out cross-validated AUC of a logistic model on X."""
    X = np.asarray(X, dtype=float).reshape(len(y), -1)
    preds = np.empty(len(y))
    model = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=5000, C=1.0))
    for train, test in LeaveOneOut().split(X):
        if len(np.unique(y[train])) < 2:
            preds[test] = np.nan
            continue
        model.fit(X[train], y[train])
        preds[test] = model.predict_proba(X[test])[:, 1]
    ok = ~np.isnan(preds)
    return auc_and_test(preds[ok], y[ok])[0]


def describe_outcome(df, y, name):
    print(f'\n{"=" * 72}\nOUTCOME: {name}   ({y.sum()}/{len(y)} positive)\n{"=" * 72}')

    for col, label in PREDICTORS.items():
        x = df[col].astype(float)
        auc, p = auc_and_test(x, y)
        lo, hi = auc_ci(x, y)
        print(f'\n{label}')
        if df[col].dtype == bool:
            tab = pd.crosstab(df[col], y)
            odds, fisher_p = stats.fisher_exact(tab.values)
            for lvl in (False, True):
                n = int(tab.loc[lvl].sum()); k = int(tab.loc[lvl, True])
                print(f'  {str(lvl):>5}: {k}/{n} = {k / n:5.1%}')
            print(f'  odds ratio {odds:.2f}   Fisher exact p = {fisher_p:.3f}')
        else:
            for lvl, sub in ((True, x[y]), (False, x[~y])):
                print(f'  {"positive" if lvl else "negative"}: '
                      f'median {sub.median():+.4f}  '
                      f'IQR [{sub.quantile(.25):+.4f}, {sub.quantile(.75):+.4f}]')
            print(f'  Mann-Whitney p = {p:.4f}')
        direction = ('higher' if auc > 0.5 else 'lower') + ' values -> positive'
        print(f'  AUC = {auc:.3f}  95% CI [{lo:.3f}, {hi:.3f}]'
              f'   (0.5 = no discrimination; {direction})')

    # LOO-CV AUC is not meaningful for a lone binary predictor (the model
    # emits two tied values and the sign can flip under resampling), so the
    # Fisher test above is the measure for that one.
    print('\nLeave-one-out cross-validated AUC:')
    for col, label in PREDICTORS.items():
        if df[col].dtype == bool:
            print(f'  {label:<38}     -- (binary; see Fisher exact above)')
        else:
            print(f'  {label:<38} {loo_auc(df[[col]], y):.3f}')
    combo = loo_auc(df[list(PREDICTORS)], y)
    print(f'  {"all three combined":<38} {combo:.3f}')
    pair = loo_auc(df[['disorder_probability', 'oqmd_stability']], y)
    print(f'  {"disorder + hull distance":<38} {pair:.3f}')


def threshold_table(df, y, name):
    """Success rate under simple screening rules."""
    print(f'\n{"=" * 72}\nSCREENING RULES for {name}\n{"=" * 72}')
    rules = {
        'no filter':                        pd.Series(True, index=df.index),
        'P_disorder < 0.25':                df.disorder_probability < 0.25,
        'hull distance < -0.015 eV/atom':   df.oqmd_stability < -0.015,
        'dynamically stable ground state':  ~df.dyn_unstable,
        'P_disorder < 0.25 AND dyn. stable':
            (df.disorder_probability < 0.25) & ~df.dyn_unstable,
        'P_disorder < 0.25 AND hull < -0.015':
            (df.disorder_probability < 0.25) & (df.oqmd_stability < -0.015),
        'all three':
            (df.disorder_probability < 0.25) & (df.oqmd_stability < -0.015)
            & ~df.dyn_unstable,
    }
    print(f'{"rule":<40} {"kept":>6} {"hits":>6} {"precision":>10} {"recall":>8}')
    total = int(y.sum())
    for label, mask in rules.items():
        kept = int(mask.sum()); hits = int((mask & y).sum())
        prec = hits / kept if kept else float('nan')
        print(f'{label:<40} {kept:>6} {hits:>6} {prec:>9.1%} '
              f'{hits / total:>7.1%}')


def loo_scores(X, y):
    """Leave-one-out cross-validated success probability for each target."""
    X = np.asarray(X, dtype=float).reshape(len(y), -1)
    out = np.empty(len(y))
    model = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=5000, C=1.0))
    for train, test in LeaveOneOut().split(X):
        model.fit(X[train], y[train])
        out[test] = model.predict_proba(X[test])[:, 1]
    return out


def precision_at_k(scores, y, ks):
    """Hit rate among the top-k ranked candidates (ties broken by index)."""
    order = np.argsort(-scores, kind='stable')
    hits = np.cumsum(y[order])
    return {k: (hits[k - 1], hits[k - 1] / k) for k in ks}


def average_precision(scores, y):
    """Mean of precision@k evaluated at every rank holding a true positive."""
    order = np.argsort(-scores, kind='stable')
    yy = y[order]
    hits = np.cumsum(yy)
    prec = hits / np.arange(1, len(yy) + 1)
    return prec[yy].mean()


def budget_analysis(df, y):
    """Threshold-free view: rank candidates, then spend a fixed budget."""
    print(f'\n{"=" * 72}\nFIXED-BUDGET RANKING (no thresholds)\n{"=" * 72}')
    print('Leave-one-out cross-validated ranking; "successes in top k melts".')

    models = {
        'hull distance alone':        ['oqmd_stability'],
        'disorder alone':             ['disorder_probability'],
        'disorder + hull':            ['disorder_probability', 'oqmd_stability'],
        'disorder + hull + dynamic':  ['disorder_probability', 'oqmd_stability',
                                       'dyn_unstable'],
    }
    ks = [5, 8, 10, 15, 20]
    base = y.mean()

    header = f'{"ranking":<28}' + ''.join(f'{f"top {k}":>10}' for k in ks) + f'{"AP":>8}'
    print('\n' + header)
    print('-' * len(header))
    for label, cols in models.items():
        sc = loo_scores(df[cols], y)
        pk = precision_at_k(sc, y, ks)
        row = f'{label:<28}'
        for k in ks:
            hits, prec = pk[k]
            row += f'{f"{int(hits)}/{k}":>10}'
        row += f'{average_precision(sc, y):>8.3f}'
        print(row)

    # A single raw feature needs no model at all - just sort by it.
    sc = -df['oqmd_stability'].to_numpy()
    pk = precision_at_k(sc, y, ks)
    row = f'{"sort by hull depth (no fit)":<28}'
    for k in ks:
        row += f'{f"{int(pk[k][0])}/{k}":>10}'
    row += f'{average_precision(sc, y):>8.3f}'
    print(row)

    row = f'{"random order (expected)":<28}'
    for k in ks:
        row += f'{f"{base * k:.1f}/{k}":>10}'
    print(row + f'{base:>8.3f}')
    print(f'\nBase rate = {base:.1%}.  AP = average precision '
          f'(threshold-free; {base:.3f} = no skill).')


def main():
    df = load_classified_samples()
    df['dyn_unstable'] = df['formula'].map(load_dynamic_stability())
    df = df.reset_index(drop=True)

    print(f'{len(df)} classified targets  '
          f'(SS {(df["XRD Result"] == "SS").sum()}, '
          f'MP {(df["XRD Result"] == "MP").sum()}, '
          f'P {(df["XRD Result"] == "P").sum()})')
    print('\nOutcome vs dynamic stability of the ground-state polymorph:')
    print(pd.crosstab(df['XRD Result'], df['dyn_unstable'],
                      rownames=['XRD'], colnames=['dyn. unstable'], margins=True))

    success = (df['XRD Result'] == 'P').to_numpy()
    describe_outcome(df, success, 'predicted structure found (P)')
    threshold_table(df, pd.Series(success, index=df.index),
                    'predicted structure found (P)')

    budget_analysis(df, success)

    disordered = (df['XRD Result'] == 'SS').to_numpy()
    describe_outcome(df, disordered, 'formed a solid solution (SS)')


if __name__ == '__main__':
    main()
