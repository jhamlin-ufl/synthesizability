# Candidate lists

The 34 targets came from three candidate lists (paper, Methods, "Candidate selection").
The lists are frozen: they define the experiment and are not rerun, because the
databases have changed since the targets were chosen and a rerun would pick different
targets.

| File | Contents |
|---|---|
| `targets.csv` | the 34 targets: list, formula as printed in the paper, sample number, and for List C the entry in the source list |
| `list_C_arcmelt.csv` | List C as selected (15 entries); five of them were made |
| `list_C_high_tc.csv` | the parallel List C selection with a predicted Tc above 2 K (56 entries); three of them were made |

The sample records (`samples/*/buttons/*/STATUS`) use older names for the lists: "B (two active
binaries, no real ternary phase)" is List A, "C (two active binaries, single real ternary
phase)" is List B, and "James's (thermodynamically stable, low-disorder superconductor
candidates from the Alexandria database)" is List C.

## List A (10 targets)

Criteria (Methods): ternary entries of MP and OQMD within 100 meV/atom of the hull, in
systems with reported binary compounds on two of the three binary edges and no reported
ternary compound in the Materials Platform for Data Science (MPDS).

Not traced: the full candidate list, the code or spreadsheet that produced it, and
the selection date. Only the 10 targets that were made are known (`targets.csv`, from the
`List:` line of each sample's STATUS file). The MPDS criterion cannot be rerun without an
MPDS licence.

## List B (16 targets)

Criteria (Methods): entries that OQMD listed on the hull at the time of selection, in
systems with binary compounds on all three edges and exactly one reported ternary
compound.

Not traced: the full candidate list, the code that produced it, and the selection
date. Only the 16 targets that were made are known. The earliest List B scans in
`samples/` (0486, 0487) are dated 2025-09-02/03, so the list predates them. Sample 0493
(TiPt5Si) has no `List:` line in its STATUS file; it is assigned to List B as in the
paper's Table 1 and with its sister sample 0494.

## List C (8 targets)

Criteria (Methods): ternary compounds on the PBEsol hull of the Alexandria subset
computed by Gao et al. (Nat. Commun. 16, 8253 (2025)), with P_disorder < 0.2
(Jakob et al. classifier); compositions at which OQMD has an entry computed from an
ICSD record removed; compounds with Tc, Pm, Th, Pa, U, Np, Pu, Tl or Ac removed; the
magnetic rare earths Ce to Yb removed; elements on hand and suitable for arc melting.
LaRuC2, NbNi2B2 and ScCrC2 came from the parallel selection that also required a
predicted superconducting Tc above 2 K (`list_C_high_tc.csv`).

Source: the files `02_arcmelt_candidates.csv` and `01_high_tc_candidates.csv` of the
group's `disorder-classification` repository (script
`scripts/complete_alexandria_pipeline_with_cached_disorder.py`), copied unchanged and
renamed. Their file dates are 2025-10-17; the first List C scan (0503) is dated
2025-11-13. All eight List C targets are entries of these two files (`targets.csv`
gives the Alexandria ids). Columns: `mat_id` Alexandria id, `formula` as stored,
`tc_010` Allen-Dynes Tc (K) from Gao et al., `lambda` and `wlog` (K) electron-phonon
coupling and log-average phonon frequency, `e_above_hull` and `e_form` (eV/atom, PBEsol),
`disorder_probability` P_disorder, `total_mag` (mu_B), `reduced_formula`.

Not traced: the selection date beyond the file dates above. The pipeline needs the
Alexandria 2023.12 PBEsol hull, the Gao et al. phonon data and a local OQMD database for
the ICSD filter; it is kept in that repository, not here.
