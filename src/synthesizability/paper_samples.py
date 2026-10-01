"""The 33 arc-melted targets that the paper reports on, and their XRD outcomes.

Every paper figure that talks about "the 33 samples" should get its sample list
from here, so that a reclassification or an added caveat lands in every figure at
once rather than in whichever script happens to be edited.

The set is defined as: the samples in the Final Analysis sheet that are not from
the Diffusion Model list (those are the BCC superconductor candidates reported
separately) and that carry an XRD outcome call.  The outcome itself comes from
the Sheet1 'XRD Result' column, overridden by ``RECLASSIFY`` below.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

SYNTH_CSV = Path('data/processed/synthesis_data.csv')
FINAL_CSV = Path('data/temp/Arc Melted Samples - Final Analysis.csv')
SHEET1_CSV = Path('data/temp/Arc Melted Samples - Sheet1.csv')

# XRD outcome codes, in the order the paper discusses them: the disordered
# solid solutions first, then the multiphase mixtures, then the samples that
# did form something like the predicted structure.
OUTCOMES = ('SS', 'MP', 'P')

OUTCOME_LABEL = {
    'SS': 'Solid Solution',
    'MP': 'Multi-Phase',
    'P': 'Predicted Structure Found',
}

# Outcome overrides applied on top of the Sheet1 'XRD Result' column, keyed by
# sample number.  Sheet1 stays the record of the original Jade call; these are
# the reclassifications agreed after re-reading the whole-pattern fits.
RECLASSIFY = {
    # 0538 CoSiY: the 88.8 wt% phase is the Pnma (oP12) polymorph, not the
    # predicted P4/nmm ground state (OQMD 1368530, -42.9 meV).  The nearest
    # OQMD Co1Si1Y1 polymorphs sit at +31 meV, ~74 meV above the prediction,
    # and the fit also carries YCo2Si2 (8.7%) and Si (2.4%).
    538: 'MP',
    # 0472 ScCu4Pd: inconsistent patterns across three syntheses and two
    # instruments, no phase identified; not a solid solution (review 2026-09).
    472: 'MP',
    # 0483 La5Al3Co2: formed La4Co3Al3 (Bragg-R 7%) plus, by mass balance,
    # La3Al and LaAl; multiphase (review 2026-09-29).
    483: 'MP',
    # 0556 Gd3CoSn6: the phase that formed is the known CeNiSi2-type
    # GdCo0.27Sn2, of which the predicted Gd3CoSn6 is an ordered 1x1x3 child
    # (Co on one third of the partially occupied sites).  A disordered variant
    # of the prediction, i.e. a solid solution (review 2026-09-30).
    556: 'SS',
    # 0512 C2CrSc: no call in Sheet1; the scan is sharp bcc Cr plus broad
    # hcp-Sc lines, no sign of the predicted ScCrC2 (which is itself an ICSD
    # compound, OQMD 14516).  Multiphase (review 2026-10-01).
    512: 'MP',
}

# Short qualifiers printed next to the formula in figures.
CAVEAT = {
    # 0507 CuSnY: P6_3mc match is excellent (91.2 wt%, Bragg-R 4.63%) but the
    # phase is PDF 04-014-1199, Sebastian et al. (2006) — a rediscovery.
    507: 'previously known',
    # 0539 CoSi3Y: ~80 wt% of a YCoSi3 polymorph 8 meV above the predicted
    # I4mm ground state (OQMD 1365767) after annealing; a partial success.
    539: 'partial: +8 meV polymorph',
    # 0556 Gd3CoSn6: known disordered parent of the predicted structure.
    556: 'disordered variant, known',
    # 0538 CoSiY: known Pnma polymorph, not the predicted one.
    538: 'known polymorph',
}

# Samples drawn with a half-filled marker in the disorder-vs-stability figure:
# the predicted structure type formed, but as a polymorph that is not the most
# stable one, or not as the majority phase.  They keep their outcome code.
PARTIAL = {539}

# Descriptive XRD outcome per sample for the sample table: what formed, in
# words, or why no call could be made.  Each entry must fit one table line
# (about 70 characters at footnotesize): RevTeX's tabular cannot wrap text in a
# p{} column on current TeX Live.
TABLE_COMMENT = {
    465: 'B2-like; cut-face scan is textured, powder scan pending',
    466: 'Fully disordered bcc (A2); unchanged by annealing',
    469: 'Not the predicted phase; unidentified',
    472: 'Not the predicted phase; unidentified, varies between syntheses',
    473: 'Disordered B2: Sc on one site (0.94), Cu/Pd on the other; annealed',
    474: 'Disordered B2: Sc on one site (0.89), Cu/Pd on the other',
    477: 'fcc solid solution; remake pending (16\\% mass loss)',
    478: 'Disordered L1$_2$: Sc on the corner site (0.77), Pd/Au on the faces',
    479: 'MoSi$_2$-type Au$_2$(Al,Mn), Al/Mn mixed; 19$^\\circ$ line absent',
    481: 'fcc (Ni,Pt,Ge) + Pt germanide; 16$^\\circ$ line absent',
    535: 'Disordered B2: Al on one site (0.92), Fe/Co on the other',
    545: 'Disordered B2: Al on one site (0.92), Fe/Co on the other',
    556: 'Known GdCo$_{0.27}$Sn$_2$, disordered parent of the prediction',
    475: 'Phases unidentified',
    476: 'Phases unidentified',
    480: 'GdPtGe + GdPt + Ge',
    483: 'La$_4$Co$_3$Al$_3$ + La--Al binaries',
    491: 'Phases unidentified',
    492: 'Phases unidentified',
    494: 'TiSi + Pt-rich fcc + Pt silicide',
    510: 'Phases unidentified',
    512: 'Cr + Sc (broad lines); carbide check pending; prediction is a known ICSD phase',
    513: 'La--Ru and La--C phases + graphite',
    514: 'Four La--Ni--Si phases',
    537: 'Gd$_3$Co$_4$Sn$_{13}$ + CoSn + GdCo$_{0.27}$Sn$_2$ + Sn',
    538: 'Known $Pnma$ polymorph of YCoSi (89\\%), 74 meV above prediction',
    540: 'Three Y--Al--Cu phases',
    553: 'GdFe$_{0.5}$Ge$_2$ + Gd$_{0.5}$Fe$_3$Ge$_3$ + Ge',
    554: 'GdFe$_{0.5}$Ge$_2$ + Gd$_{0.5}$Fe$_3$Ge$_3$ + GdFe$_2$Ge$_2$',
    555: 'GdFe$_{0.5}$Ge$_2$ + GdFe$_2$Ge$_2$ + Gd$_{0.5}$Fe$_3$Ge$_3$',
    493: 'Predicted structure + Pt-rich fcc impurity; no PDF entry',
    507: 'Predicted structure (91\%); known compound (ICSD; OQMD 26268)',
    539: 'Polymorph of the prediction, 8 meV above it (80\\%, annealed)',
    561: 'Predicted structure + impurities; no PDF entry',
}


def load_classified_samples(synth_csv=SYNTH_CSV, final_csv=FINAL_CSV,
                            sheet1_csv=SHEET1_CSV):
    """The 33 targets with an XRD outcome, ordered SS, then MP, then P.

    Returns the ``synthesis_data`` rows with two columns added: ``XRD Result``
    (one of ``OUTCOMES``, post-``RECLASSIFY``) and ``caveat`` (from ``CAVEAT``,
    NaN where there is none).
    """
    final = pd.read_csv(final_csv)
    sheet1 = pd.read_csv(sheet1_csv)
    synth = pd.read_csv(synth_csv)

    final['Sample Number'] = pd.to_numeric(final['Sample Number'], errors='coerce')
    sheet1['Sample Number'] = pd.to_numeric(sheet1['Sample Number'], errors='coerce')
    sheet1 = sheet1.dropna(subset=['Sample Number'])
    sheet1['Sample Number'] = sheet1['Sample Number'].astype(int)

    nums = final['Sample Number'].dropna().astype(int).tolist()
    sub = synth[synth['sample_number'].isin(nums) &
                (synth['prediction_list'] != 'Diffusion Model')].copy()
    sub = sub.join(sheet1.set_index('Sample Number')[['XRD Result']],
                   on='sample_number')

    unknown = set(RECLASSIFY) - set(sub['sample_number'])
    if unknown:
        raise ValueError(f'RECLASSIFY names samples not in the figure set: '
                         f'{sorted(unknown)}')
    sub['XRD Result'] = sub['XRD Result'].mask(
        sub['sample_number'].isin(RECLASSIFY),
        sub['sample_number'].map(RECLASSIFY))

    sub = sub[sub['XRD Result'].isin(OUTCOMES)].copy()

    unknown = set(CAVEAT) - set(sub['sample_number'])
    if unknown:
        raise ValueError(f'CAVEAT names samples not in the figure set: '
                         f'{sorted(unknown)}')
    sub['caveat'] = sub['sample_number'].map(CAVEAT)

    sub['_order'] = sub['XRD Result'].map({o: i for i, o in enumerate(OUTCOMES)})
    return (sub.sort_values(['_order', 'sample_number'])
               .drop(columns='_order')
               .reset_index(drop=True))
