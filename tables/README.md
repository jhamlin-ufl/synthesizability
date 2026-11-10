# Tables

Paths in these tables are relative to the repository root.

## `paper_samples.csv`: one row per target, in the order of Table 1

| Column | Meaning |
|---|---|
| `table_row` | row in Table 1 |
| `formula` | formula as printed in the paper; also the folder name in `samples/` |
| `list` | candidate list A, B or C (`lists/README.md`) |
| `category`, `category_name` | XRD outcome: 1 On-hull structure; 2 Above-hull polymorph; 3 Disordered parent of on-hull structure; 4 Other solid solution; 5 Multiphase, on-hull structure not formed; 6 On-hull structure did not form but structure(s) unresolved |
| `alt_category` | a second category the data cannot exclude (split marker in Fig. 2) |
| `confidence` | firm, probable or tentative |
| `reported` | `yes` structure reported before at this composition; `composition` composition reported without a structure; `series` ordered member of a reported series whose known member is partly disordered; `no` |
| `sample` | sample number of the target |
| `record_sample` | button whose scan gives the outcome |
| `record_scan` | the scan used |
| `annealed` | `yes` if the scan used is of a sample annealed at 1000 C for one week |
| `instrument`, `mount` | instrument and mount of the scan used |
| `mass_loss_pct` | arc-melting mass loss of `record_sample` (%), from its SYNTHESIS file |
| `other_samples` | the other buttons of the target and their role |

## `outcomes.csv`: one row per target

| Column | Meaning |
|---|---|
| `sample` | sample number of the target |
| `one_or_three` | `yes` where category 1 cannot be excluded for a category 3 target |
| `outcome`, `outcome_tex` | the outcome text of Table 1 (plain and LaTeX) |
| `phases` | phases and weight fractions of the main fit (Table 1 phase table) |

## `buttons.csv`: one row per arc-melted button

| Column | Meaning |
|---|---|
| `sample`, `folder` | button number and its folder |
| `formula`, `target_sample` | the target the button belongs to |
| `role` | `target`; `earlier attempt` (replaced by the target button); `remake` (made after it); `second button` |
| `mass_loss_pct` | as written in SYNTHESIS |
| `n_scans` | number of scans of the button |
| `note` | why the button was made or replaced. "Weighed composition off target" gives the largest deviation of the weighed masses from the target stoichiometry |
| `in_table_s1` | whether Table S1 lists the button |

## `scans.csv`: one row per powder pattern

| Column | Meaning |
|---|---|
| `sample`, `file` | button number and scan file |
| `instrument` | Siemens D500 (`.txt`, Cu K-alpha, 2/2/0.05 deg slits) or PANalytical X'Pert MPD of the UF Nanoscale Research Facility (`.xy`) |
| `mount` | see below |
| `state` | as cast or annealed (1000 C, one week); 0615 also has a scan after a further 500 C, 12 h anneal |
| `two_theta_start`, `two_theta_stop`, `step` | scan range and step (deg) |
| `date` | measurement date from the file header (blank for `.xy`) |
| `scan_of_record` | `yes` for the 34 scans used |
| `note` | caveats |

Mounts: scans with `_ZB` in the name are powder on a zero-background Si plate. All PANalytical
scans are of bulk pieces (0474, 0478 and 0465 cut faces, 0477 an uncut face). 0492 was scanned as
a bulk chunk on the Siemens. For the other Siemens scans the mount was not recorded. They are
taken to be powder fixed with grease on a PMMA slide, the other mount given in the Methods.

## `si_fits.csv`: one row per fit shown in the Supplement

| Column | Meaning |
|---|---|
| `table_row`, `formula` | Table 1 row and formula |
| `sample` | button whose scan was fitted |
| `role` | `FoR` main fit; `C1`, `C2` contrast models; `ascast` the as-cast scan of an annealed target |
| `label` | model name in the figure legend |
| `fit_dir` | the fit folder |
| `target_phase` | phase drawn in green (a key of `ticks.json`); blank if no phase is at the target composition |
| `zoom` | zoom windows of panel (c), `lo-hi;lo-hi` in deg |
| `sticks` | a CIF whose stick pattern replaces the model curve in the zoom panels (none at present) |
| `claim` | what the fit shows |

## `lists/`

The three candidate lists the targets came from, with their criteria (`lists/README.md`).
