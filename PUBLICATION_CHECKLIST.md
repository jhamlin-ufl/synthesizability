# Pre-publication checklist

Started 2026-10-01. Things that have to be tracked down, decided or re-done
before the synthesizability paper is submitted. Tracked in git on purpose so it
is not lost; tick items off here rather than elsewhere. Items are grouped,
roughly in dependency order. Related working notes (XRD review, student action
items, occupancy analysis) are in `data/temp/`, which is gitignored.

## State at the end of the 2026-10-01 session (read first)

- **Committed on branch `paper-figures-tables`** (2026-10-01, commit
  6935d39, 684 files; not pushed, not merged). Left uncommitted on purpose:
  `reference_papers/` (publisher PDFs), `notes_on_xrd.txt` (candid voice
  notes) and `synthesizability_study_structures.zip` (origin unclear). Still
  to do: tag the database state (`git tag db-snapshot-A`), push, and open the
  PR when the figures settle.
- **Overleaf** (`~/ufdb/Apps/Overleaf/Synthesizability Project/`): main.tex,
  supplementary.tex, `figures/{xrd_representative,xrd_appendix_SS,MP,P,
  disorder_vs_stability}.pdf`, `tables/{sample_table,database_table}.tex`
  are copies of `results/publication_ready/` as of 2026-10-01 ~12:30.
  `hull_cross_database.pdf` is NOT yet in the paper. Overleaf receives edits
  only while the UFL Dropbox client is running (it died after a self-update
  on 2026-10-01 morning; the user restarts it with their own procedure).
- **Working notes (gitignored, `data/temp/`)**: `xrd_notes_consolidated.md`
  (per-sample XRD conclusions, authoritative), `student_action_items_
  2026-09-30.md` (what the student was asked to fit; item 6 = remake 0512),
  `jade_cifs/` + `jade_cifs_2026-09-30.zip` (CIF sets given to the student),
  `known_structure_audit_2026-10-01.csv`, occupancy-analysis scripts/PDFs.
- **Next session: reviewing the student's Jade fits.** For each returned fit:
  the exports belong in `data/raw/<sample>/` as `<nnnn>_XRD_fit.txt`
  (Angle, I(o), I(c), I(d), I(b)) and `<id>.wpf.txt`; check that the phases
  are restricted to the right chemistry, that the model is the one asked for
  (table in the student file), compare refined site occupancies and lattice
  constants with ours (`SOLUTION` in `src/synthesizability/xrd_figure.py`:
  466 A2 a=3.258; 473 B2 Sc 0.94 a=3.277; 474 B2 Sc 0.89 a=3.271; 477 fcc
  3.968; 478 L1_2 corner 0.77 a=3.980; 479 C11b 3.333/8.897; 481 fcc 3.587;
  535/545 B2 Al 0.92 a=2.864/2.868) and the Bragg-R values. To adopt a fit:
  point the sample's `SOLUTION` entry at `Jade()` (or keep the model and
  update its numbers), update `TABLE_COMMENT` in `paper_samples.py`, run
  `snakemake results/publication_ready/{xrd_representative.pdf,
  xrd_appendix_SS.pdf,sample_table.tex}` and copy to Overleaf.
  Anneal/remake decisions still open: 0477 ScPd5Au2 remake (16 % loss),
  0512 C2CrSc remake, Cu5GePd6 powder scan, 0480 GdGePt2 Siemens scan push.
- **Dashboard** has not been rebuilt since the Alexandria switch
  (`snakemake results/dashboard/index.html`; then the usual refresh PR).

## Databases

- [ ] **Snapshot A is recorded** in `data/external/snapshot_provenance.json`
      (OQMD v1.7 local dump, 1.32M entries; MP API April 2026; Alexandria
      OPTIMADE Sep 2026 for PBEsol; Alexandria PBE from the bulk 2025.07.02
      release, 58 chunk files, listing hash recorded; re-recorded 2026-10-01 12:00). When this state is committed, tag it
      (`git tag db-snapshot-A`) so it can be checked out later.
- [ ] Re-pull all three databases shortly before submission into a parallel
      cache ("snapshot B"), without overwriting A: OQMD (newer dump or API),
      MP API, Alexandria OPTIMADE and the newest convex-hull release.
- [ ] Diff A against B per target: hull distance, on-hull status, lowest-entry
      structure, existence of the predicted entry. This is the "databases are
      in flux" paragraph; Alexandria already moved four targets off its hull
      between 2023.12 and 2025.07 and has Hf2MoIr on opposite sides of the hull
      in the 2025.07 file and the current API.
- [ ] Decide which snapshot the paper's numbers quote (recommendation: B, with
      A to B changes reported) and state the versions and dates in Methods.
- [ ] Methods sentence: all hull distances are recomputed from the cached
      entries by removing every entry at the target composition and building
      the hull of the rest (`compute_hull_distances.py`); OQMD's stored
      stability column is not used (Hf2MoIr: stored -18.6 vs recomputed
      -12.0 meV/atom; the live OQMD API reports 0 for all hull phases).
- [x] Alexandria structures for compositions not on its hull: obsolete since
      2026-10-01, the bulk-release extractor writes every entry at the target
      compositions (245 structures, on and off the hull).
- [ ] YCuSn: the OQMD "prediction" (entry 1468383) is a prototype duplicate of
      the ICSD-derived OQMD entry 26268 (same structure, energies equal to
      0.1 meV/atom). Say so in the text; it is a database-rediscovery case,
      not a prediction. Consider an `icsd` screen in the target-selection code
      for future lists.
- [ ] CoSiY: OQMD and Alexandria predict P4/nmm; MP's ground state is the
      Pnma phase (mp-1207632) that actually formed; OQMD has Pnma at
      +31 meV. Worth a sentence.
- [ ] Sc2Pd5Au: OQMD and Alexandria Pbam, MP Amm2, experiment disordered
      L1_2. Worth a sentence.
- [ ] **Alexandria source decision.** The repo's Alexandria hull distances
      (`data/external/alexandria_*_ternary_phases`, used by
      `compute_hull_distances.py` and shown in Table S2) come from an OPTIMADE
      API pull of 2026-04-08, repaired in place on 2026-09-28. That pull was a
      false start: the API was unreliable and the pull is incomplete (see next
      item). The candidate-selection work (`~/repos/candidate_materials`,
      2026-09-28) switched to the downloaded convex-hull release files in
      `~/data/alexandria/` (PBE 2025.07.02; PBE and PBEsol 2023.12.29), which are
      versioned, checksummed in `snapshot_provenance.json`, and offline; this
      repo now reads them only for structures
      (`extract_alexandria_hull_structures.py`). Plan: rebuild the Alexandria
      hull column from the 2025.07 release file. A release file holds only
      on-hull entries, so it gives the hull and the structure of every stable
      phase but not the energy of an off-hull target; for the off-hull targets
      take the target composition's own entries from a per-composition OPTIMADE
      query or from the full Alexandria download, and record which in Methods.
      Until this is done the Alexandria column is the incomplete OPTIMADE pull.
- [ ] Alexandria OPTIMADE cache vs 2025.07 hull file disagree on coverage:
      TiGePt5, TiSiPt5, ZrSiPt5, Al4FeCo3 and Sc2CuPd are on the 2025.07 hull
      with structures, but the OPTIMADE pull has no entry (blank hull distance)
      or a different verdict for some of them. Re-check the OPTIMADE query
      (composition filter) at the refresh; until then Table S2 shows blanks
      next to "yes" for those rows.
- [ ] Decide whether Alexandria PBEsol appears anywhere (recommendation: no;
      different functional, 14/33 coverage).

## Samples and classification

- [x] 512 C2CrSc: in, as multiphase (decided 2026-10-01), making the set 34.
      Sheet1 still has no XRD call for it; `RECLASSIFY` supplies MP. Note the
      predicted ScCrC2 (OQMD 14516) is ICSD-derived: a second rediscovery
      target after YCuSn, and one that did not form. Student refit pending
      (Cr + Sc + carbides).
- [ ] 568 Cu5GePd6 and 566 ScCu4Pd (remakes) appeared in the old table in place
      of 465 and 472. Paper set follows the Final Analysis sheet; confirm.
- [ ] Marker for 507 (predicted phase, previously known) and 538 (known
      polymorph) in the disorder-vs-stability figure: half-filled or not.
      `PARTIAL` in `paper_samples.py` currently holds only 539.
- [ ] Move the two sample sheets out of gitignored `data/temp` so the figure
      and table rules can be pushed; until then `snakemake` fails for anyone
      else.
- [ ] Discussion text still quotes the old counts ("27/28 solid solutions",
      "six(?) samples"); rewrite against the table (12 SS, 18 MP, 4 P; 34 targets).

## XRD data and fits still outstanding (see `data/temp/student_action_items_2026-09-30.md`)

- [ ] 0480 GdGePt2: Siemens 5-120 deg scan used for the Jade fit is not in the
      repo; dashboard shows only the Panalytical scan.
- [ ] 0466 Hf2MoIr annealed Jade fit files (HfMo3 + Mo) not in the repo.
- [ ] Refits with the refined-occupancy CIFs (B2/L1_2/A2/C11b) and the
      multiphase refits (La5Al3Co2 + La3Al; CoSi3Y annealed without the "Si12"
      phase; TiSiPt5 / ZrSiPt5 with the predicted structure + Pt-rich fcc;
      C2LaRu without the cubic carbon; CoSiY improved; Ti2SiPt8 without Zr).
      Jade exports (JPG, wpf, wrk) pushed per sample.
- [ ] Cu5GePd6: powder scan (cut-face scan is (110)-textured; intensities
      unusable). Then B2 vs predicted Cmmm fit.
- [ ] ScPd5Au2: remake (16 % mass loss) and 5 deg scan.
- [ ] Sc2Pd5Au: 5 deg scan to reach the 15.7 and 17.5 deg predicted lines.
- [ ] Grinding test on one of ScCuPd2 / ScCu4Pd / Sc3CuPd2 (Panalytical bulk vs
      Siemens powder scans disagree).
- [x] Alexandria PBE source switched (2026-10-01): the per-space PBE cache is
      now built from the bulk 2025.07.02 release (58 files, 3.4 GB, in
      ~/data/alexandria/pbe_2025.07.02; `extract_alexandria_release_phases.py`),
      not from OPTIMADE. Reason: the OPTIMADE PBE server serves only part of
      the database (missed on-hull TiPt5Si, TiPt5Ge; id ranges stop early in
      several spaces) and was returning HTTP 500 on 2026-10-01. PBEsol still
      comes from OPTIMADE (live counts match the cache exactly). The old
      OPTIMADE PBE cache is not kept in the repo; a copy sits in the session
      scratchpad only. Consequences: Alexandria hull distances changed for
      about a third of the targets (e.g. Hf2MoIr -2.9 -> +3.7, ScCu4Pd -18.0 ->
      +3.9, NiGe2Pt 236 -> 17.9 meV), Alexandria structure agreement is now
      assessed for 30/34 targets (21 same), Table S1 rows for Ti2Pt8Ge, TiPt5Ge
      and TiPt5Si are complete.
- [ ] Alexandria validation: our recomputed Alexandria hull distances agree
      with Alexandria's own `e_phase_separation` to < 0.5 meV for 25 of 30
      targets; five on-hull targets differ by 7-17 meV with ours less
      negative (TiPt5Ge, TiPt5Si, ZrPt5Si, YCoSi, GdFeGe3). Most likely
      Alexandria's stored value predates newer competing entries; ours is
      self-consistent with the release. State which is used (ours) in the
      Methods, or ask the Alexandria group.
- [ ] Dashboard needs a rebuild against the new Alexandria PBE cache (schema
      unchanged, extra fields added; per-space CIFs now exist for PBE).
- [ ] C2CrSc (0512): Pöttgen et al., J. Solid State Chem. 119, 324 (1995)
      (now in `reference_papers/Pottgen1995_ScCrC2_JSSC119_324.pdf`) shows the
      compound forms single-phase directly from the arc melt when the charge
      is a pressed pellet of filings/powder/graphite flakes; annealing only
      orders it (alpha form). Our button (solid pieces, low current, short
      melts) is Cr + Sc with <= 5 wt% carbide. Decide remake vs anneal (student
      item 6); then settle the Table I comment ("carbide check pending") and
      the Discussion wording (known compound missed by our standard recipe).
      Sc-Cr-C subsolidus diagram: Artyukh et al., Poroshkovaya Metallurgiya
      (1997), not yet obtained.
- [ ] ICSD leak audit (2026-10-01, from the OQMD cache `icsd` flag at each
      target composition): three of the 34 targets have ICSD-derived OQMD
      entries at the composition. C2CrSc: the selected entry 14516 itself is
      ICSD-derived (list J, built from Gao 2025 / Alexandria, where no OQMD
      ICSD check was applied). CuSnY: the selected entry 1468383 is a
      non-ICSD duplicate of ICSD entry 26268 (list J). GdFeGe2: ICSD entry
      30184 (Cmcm, 8 atoms, +59 meV in OQMD) is the CeNiSi2-type defect phase
      that dominates the 0555 fit (GdFe0.52Ge2, 61 wt%); the selected P2_1/m
      entry 1377698 is a different structure (list C, screened by literature
      search). Decide how to present these in Methods/Discussion; the
      Methods sentence "not already present in the ICSD" for list J needs
      qualifying. Verified 2026-10-01 against the PDF card: the predicted P2_1/m cell is
      the fully occupied CeNiSi2-type framework (homogeneity range of the
      real phase 0.25 <= x <= 0.46); 0555 formed the equilibrium
      GdFe_xGe2 + GdFe2Ge2 (+ GdFe6Ge6) assemblage of the 800 C section.
      Sources in reference_papers/ (Jemmali 2010, 2016).
- [ ] OQMD provenance flags, what they mean for the three "known" cases
      (2026-10-01): OQMD's `icsd` keyword marks entries whose input structure
      was an ordered ICSD record (OQMD imports no partial occupancies), not
      that a stoichiometric compound exists. ScCrC2 14515/14516: genuine
      single-crystal ICSD records of the ordered compound. CuSnY 26268:
      genuine ICSD record; our selected 1468383 is a prototype-library
      (LiGaGe_type) entry that OQMD itself marks `hold: duplicate`. GdFeGe2
      30184: ICSD record 631913 is an idealized fully occupied cell of the
      defect phase GdFe_xGe2 (x <= 0.46); no stoichiometric compound exists.
      OQMD's `entries.duplicate_of_id` resolves the holds: 1468383 (CuSnY) ->
      26268 = icsd-416544 (an ICSD entry); 1991578 (ScCu4Pd, 0472) ->
      1888746, a Griesemer/Armiento hypothetical structure with no ICSD
      keyword, so 0472 remains a genuine prediction. The pipeline reads
      neither `hold` nor `duplicate_of_id`. Audit of all 34 (2026-10-01, data/temp/known_structure_audit_2026-10-01.csv,
      ad hoc; make it a Snakefile rule if the number goes in the paper):
      selected entry `icsd` or duplicate_of -> `icsd`: 2/34 known (C2CrSc,
      CuSnY), 32/34 without an ICSD structure; any ICSD entry at the
      composition: 3/34 (adds GdFeGe2's idealized defect-phase cell);
      Materials Project theoretical=False at the composition: 1/34
      (C2CrSc, mp-4992 and mp-1188534, both the predicted structure); MP
      has no entry at all for 23/34 compositions and marks its YCuSn and
      YCoSi entries theoretical. Adopted criterion (PI, 2026-10-01): a target counts as a novel
      composition when NO OQMD entry at that composition carries the `icsd`
      keyword: 31/34. The three flagged ones are stated individually:
      C2CrSc and CuSnY are known compounds at that composition; GdFeGe2 is
      flagged only via the idealized full-occupancy record of GdFe_xGe2
      (x <= 0.46), i.e. a new ordered phase, not a rediscovery. The
      duplicate_of analysis explains how CuSnY passed selection and is
      not part of the criterion. Methods should state this criterion.
- [ ] Sc3CuPd2 as-cast vs annealed B2 order: only a Rietveld refinement of both
      scans can say whether order increased; peak-vector estimates are within
      their uncertainty of each other.

## Figures and tables

- [ ] Figure 1 panel choice once refits land (TiSiPt5 gets a fit; consider
      whether ZrSiPt5 becomes the better novel example).
- [ ] Figure 1 is a float page anchored after the fourth Discussion paragraph
      so the table precedes it; move the anchor when the Results text grows.
- [ ] Figure 1 caption: predicted patterns at the PBE lattice; c/a note for
      YCuSn. Keep in sync if panel (e) changes.
- [ ] Supplement: XRD grids use the annealed scans for 466, 473, 539; the
      annealing figure (as-cast vs annealed for 466, 473, 474, 539) is not yet
      made.
- [ ] Occupancy analysis (`data/temp/occupancy_goodness.py`): decide whether
      it goes in the supplement; if so, move it into `src/` with a rule.
- [ ] Disorder-vs-stability figure reads `oqmd_stability` from the sample data
      (OQMD's stored column); switch it to `hull_distances.csv` so it agrees
      with the table (Hf2MoIr -18.6 vs -12.0).
- [ ] Formula element order (Allen electronegativity) is used in all figures
      and tables; confirm that is the convention wanted in the paper.

## Methods text to add

- [ ] Hull-distance construction (above) and database snapshot dates.
- [ ] Peak-vector / peak-fit occupancy method, if used.
- [ ] Annealing conditions (1000 C, one week, sealed under Ar) and the
      selection rule for annealed samples (solid-solution candidates with a
      diagnostic low-angle line and material remaining).
- [ ] Whole-pattern fitting: Jade version, which scan per sample (instrument,
      bulk vs powder), background and profile settings.
