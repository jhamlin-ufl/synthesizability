# Snakefile for synthesizability pipeline
import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Data file lists
# ---------------------------------------------------------------------------

RAW_DATA_FILES = [str(p) for p in Path("data/raw").rglob("*") if p.is_file()]
CHI_DATA_FILES = [str(p) for p in Path("data/raw").rglob("*chiAC*.txt")]
WPF_FILES = [str(p) for p in Path("data/raw").rglob("*.wpf.txt")]
TERNARY_JSON_FILES = (
    ([str(p) for p in Path("data/external/oqmd_ternary_phases").rglob("*.json")]
     if Path("data/external/oqmd_ternary_phases").exists() else []) +
    ([str(p) for p in Path("data/external/mp_ternary_phases").rglob("*.json")]
     if Path("data/external/mp_ternary_phases").exists() else []) +
    ([str(p) for p in Path("data/external/alexandria_pbe_ternary_phases").rglob("*.json")]
     if Path("data/external/alexandria_pbe_ternary_phases").exists() else []) +
    ([str(p) for p in Path("data/external/alexandria_pbesol_ternary_phases").rglob("*.json")]
     if Path("data/external/alexandria_pbesol_ternary_phases").exists() else [])
)
GENAI_CIF_FILES = [str(p) for p in Path("data/external/genai_structures").rglob("*.cif")]
XRD_JPG_FILES = [str(p) for p in Path("data/raw").rglob("*_XRD_fit*.JPG")]
RAW_CIF_FILES = [str(p) for p in Path("data/raw").rglob("*.cif")]
STATUS_FILES = [str(p) for p in Path("data/raw").rglob("STATUS")]
SYNTHESIS_FILES = [str(p) for p in Path("data/raw").rglob("SYNTHESIS")]
XRD_FILES = (
    [str(p) for p in Path("data/raw").rglob("*.xy")] +
    [str(p) for p in Path("data/raw").rglob("*.txt") if "chiAC" not in p.name]
)
DATAFRAME_INPUT_FILES = STATUS_FILES + SYNTHESIS_FILES + XRD_FILES

REFERENCE_DATA_FILES = [
    "data/external/reference/element_prices.csv",
    "data/external/reference/element_vapor_pressures.csv",
]
DISORDER_MODEL_FILES = [
    "data/external/disorder_model/hyperopt_best_model.pt",
    "data/external/disorder_model/hyperopt_config.json",
]
SUPERCON_DATA_FILES = [
    "data/external/supercon/primary.tsv",
]

# ---------------------------------------------------------------------------
# Source file lists — targeted per rule to avoid spurious reruns
# ---------------------------------------------------------------------------

SRC_IO = (
    [str(p) for p in Path("src/synthesizability/io").rglob("*.py")] +
    [str(p) for p in Path("src/synthesizability/parsers").rglob("*.py")] +
    ["src/synthesizability/formula.py"]
)

SRC_DISORDER = (
    [str(p) for p in Path("src/synthesizability/disorder_core").rglob("*.py")] +
    ["src/synthesizability/disorder.py"]
)

SRC_SUSCEPTIBILITY = [
    "src/synthesizability/susceptibility.py",
]

SRC_OQMD = [
    "src/synthesizability/oqmd.py",
]

SRC_DASHBOARD = (
    [str(p) for p in Path("src/synthesizability/dashboard_plugins").rglob("*.py")]
)

SRC_FIGURES = [
    "src/synthesizability/paper_samples.py",
    "src/synthesizability/xrd_figure.py",
    "src/synthesizability/xrd_background.py",
    "src/synthesizability/xrd_simulator.py",
    "src/synthesizability/parsers/xrd.py",
]

# Sample classification and XRD outcome calls still live in these exported
# sheets rather than in a parsed data file.
SAMPLE_SHEETS = [
    "data/temp/Arc Melted Samples - Final Analysis.csv",
    "data/temp/Arc Melted Samples - Sheet1.csv",
]

# XRD outcome classes; one appendix XRD figure is drawn per class.
XRD_OUTCOMES = ["SS", "MP", "P"]

print(f"Tracking {len(RAW_DATA_FILES)} raw data files")
print(f"  - {len(CHI_DATA_FILES)} chi data files")
print(f"  - {len(DATAFRAME_INPUT_FILES)} dataframe input files (STATUS/SYNTHESIS/XRD)")

# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

rule all:
    input:
        "data/processed/synthesis_data.csv",
        "data/processed/synthesis_data.pkl",
        "data/processed/oqmd_hull_data.csv",
        "data/external/oqmd_structures/.extracted",
        "data/external/oqmd_ternary_phases/.extracted",
        "results/susceptibility/susceptibility_real_part.pdf",
        "results/susceptibility/susceptibility_imaginary_part.pdf",
        "results/susceptibility/hc2_with_fits.pdf",
        "results/susceptibility/hc2_fit_parameters.csv",
        "results/dashboard/index.html",
        "results/publication_ready/xrd_representative.pdf",
        expand("results/publication_ready/xrd_appendix_{outcome}.pdf", outcome=XRD_OUTCOMES),
        "results/publication_ready/disorder_vs_stability.pdf",
        "results/publication_ready/sample_table.tex",
        "results/publication_ready/database_table.tex",
        "results/publication_ready/hull_cross_database.pdf",
        "data/external/snapshot_provenance.json"


checkpoint build_dataframe_for_formulas:
    input:
        script="scripts/build_dataframe.py",
        src=SRC_IO,
        data=STATUS_FILES + SYNTHESIS_FILES,
        reference=REFERENCE_DATA_FILES,
    output:
        csv="data/processed/synthesis_data_no_disorder.csv",
        formulas="data/processed/formulas.txt",
    run:
        import shutil
        import pandas as pd

        cache_path = Path("data/processed/disorder_cache.csv")
        backup_path = Path("data/processed/disorder_cache.csv.backup")

        if cache_path.exists():
            shutil.move(str(cache_path), str(backup_path))

        shell("poetry run python {input.script}")

        shutil.move("data/processed/synthesis_data.csv", str(output.csv))

        if backup_path.exists():
            shutil.move(str(backup_path), str(cache_path))

        df = pd.read_csv(output.csv)
        formulas = sorted(df['formula'].dropna().unique().tolist())
        Path(output.formulas).write_text('\n'.join(formulas) + '\n')


rule compute_disorder_cache:
    input:
        script="scripts/compute_disorder_probabilities.py",
        src=SRC_DISORDER,
        model_files=DISORDER_MODEL_FILES,
        formulas="data/processed/formulas.txt",
    output:
        cache="data/processed/disorder_cache.csv",
    shell:
        "poetry run python {input.script}"


rule compute_supercon_cache:
    input:
        script="scripts/compute_supercon_cache.py",
        supercon=SUPERCON_DATA_FILES,
        formulas="data/processed/formulas.txt",
    output:
        marker=touch("data/processed/.supercon_cached"),
    shell:
        "poetry run python {input.script}"


rule validate_oqmd_database:
    output:
        validation_marker=touch("data/processed/.oqmd_validated"),
    log:
        "logs/validate_oqmd_database.log"
    shell:
        "poetry run python scripts/validate_oqmd_database.py > {log} 2>&1"


rule query_oqmd_hulls:
    input:
        script="scripts/query_oqmd_hulls.py",
        src=SRC_OQMD,
        validation="data/processed/.oqmd_validated",
        csv="data/processed/synthesis_data_no_disorder.csv",
    output:
        csv="data/processed/oqmd_hull_data.csv",
    log:
        "logs/query_oqmd_hulls.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


rule extract_oqmd_structures:
    input:
        script="scripts/extract_oqmd_structures.py",
        src=SRC_OQMD,
        validation="data/processed/.oqmd_validated",
        hull_data="data/processed/oqmd_hull_data.csv",
    output:
        marker=touch("data/external/oqmd_structures/.extracted"),
    log:
        "logs/extract_oqmd_structures.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


rule query_ternary_phases:
    input:
        script="scripts/query_ternary_phases.py",
        src=SRC_OQMD,
        validation="data/processed/.oqmd_validated",
        csv="data/processed/synthesis_data_no_disorder.csv",
    output:
        marker=touch("data/external/oqmd_ternary_phases/.queried"),
    log:
        "logs/query_ternary_phases.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


rule extract_ternary_cifs:
    input:
        script="scripts/extract_ternary_cifs.py",
        src=SRC_OQMD,
        queried="data/external/oqmd_ternary_phases/.queried",
    output:
        marker=touch("data/external/oqmd_ternary_phases/.extracted"),
    log:
        "logs/extract_ternary_cifs.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


rule query_mp_phases:
    input:
        script="scripts/query_mp_ternary_phases.py",
        csv="data/processed/synthesis_data_no_disorder.csv",
    output:
        marker=touch("data/external/mp_ternary_phases/.queried"),
    log:
        "logs/query_mp_phases.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


rule extract_mp_cifs:
    input:
        script="scripts/extract_mp_cifs.py",
        queried="data/external/mp_ternary_phases/.queried",
    output:
        marker=touch("data/external/mp_ternary_phases/.extracted"),
    log:
        "logs/extract_mp_cifs.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


# Alexandria PBEsol: OPTIMADE API (its content matches the live server; checked
# 2026-10-01).  Alexandria PBE comes from the bulk release instead, see
# extract_alexandria_release_phases below.
rule query_alexandria_phases:
    input:
        script="scripts/query_alexandria_phases.py",
        csv="data/processed/synthesis_data_no_disorder.csv",
    output:
        pbesol_marker=touch("data/external/alexandria_pbesol_ternary_phases/.queried"),
    log:
        "logs/query_alexandria_phases.log"
    shell:
        "poetry run python {input.script} pbesol > {log} 2>&1"


rule extract_alexandria_cifs:
    input:
        script="scripts/extract_alexandria_cifs.py",
        pbesol_queried="data/external/alexandria_pbesol_ternary_phases/.queried",
    output:
        pbesol_marker=touch("data/external/alexandria_pbesol_ternary_phases/.extracted"),
    log:
        "logs/extract_alexandria_cifs.log"
    shell:
        "poetry run python {input.script} pbesol > {log} 2>&1"


# Alexandria PBE: every entry in the target chemical spaces, with structures at
# the target compositions, from the bulk 2025.07.02 release files in
# ALEXANDRIA_DIR/pbe_2025.07.02 (58 files, ~3.4 GB, from
# https://alexandria.icams.rub.de/data/pbe/2025.07.02/).  The OPTIMADE server
# serves only part of the PBE database, so it is not used for PBE.
ALEXANDRIA_DIR = os.environ.get("ALEXANDRIA_DIR", os.path.expanduser("~/data/alexandria"))
ALEXANDRIA_PBE_CHUNKS = sorted(str(p) for p in Path(ALEXANDRIA_DIR, "pbe_2025.07.02").glob("alexandria_*.json.bz2"))

rule extract_alexandria_release_phases:
    input:
        script="scripts/extract_alexandria_release_phases.py",
        formulas="data/processed/synthesis_data_no_disorder.csv",
        chunks=ALEXANDRIA_PBE_CHUNKS,
    output:
        queried=touch("data/external/alexandria_pbe_ternary_phases/.queried"),
        extracted=touch("data/external/alexandria_pbe_ternary_phases/.extracted"),
        index="data/external/alexandria_hull_structures/index.csv",
    log:
        "logs/extract_alexandria_release_phases.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


rule compute_hull_membership:
    input:
        script="scripts/compute_hull_membership.py",
        src="src/synthesizability/dashboard_plugins/ternary_phases.py",
        formulas="data/processed/synthesis_data_no_disorder.csv",
        ternary_data=TERNARY_JSON_FILES,
    output:
        csv="data/processed/hull_membership.csv",
    shell:
        "poetry run python {input.script}"


# Hull distance of each target in every database, OQMD-style (negative = below
# the hull of the other compositions), recomputed from the cached entries.
rule compute_hull_distances:
    input:
        script="scripts/compute_hull_distances.py",
        formulas="data/processed/synthesis_data_no_disorder.csv",
        ternary_data=TERNARY_JSON_FILES,
    output:
        csv="data/processed/hull_distances.csv",
    log:
        "logs/compute_hull_distances.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


rule build_dataframe:
    input:
        script="scripts/build_dataframe.py",
        src=SRC_IO,
        data=DATAFRAME_INPUT_FILES,
        reference=REFERENCE_DATA_FILES,
        disorder_cache="data/processed/disorder_cache.csv",
        oqmd_hulls="data/processed/oqmd_hull_data.csv",
        hull_membership="data/processed/hull_membership.csv",
        remake_map="data/raw/REMAKE_MAP.csv",
    output:
        csv="data/processed/synthesis_data.csv",
        pkl="data/processed/synthesis_data.pkl",
    shell:
        "poetry run python {input.script}"


rule analyze_susceptibility:
    input:
        script="scripts/analyze_susceptibility.py",
        src=SRC_IO + SRC_SUSCEPTIBILITY,
        data=CHI_DATA_FILES,
    output:
        real="results/susceptibility/susceptibility_real_part.pdf",
        imag="results/susceptibility/susceptibility_imaginary_part.pdf",
        hc2="results/susceptibility/hc2_with_fits.pdf",
        params="results/susceptibility/hc2_fit_parameters.csv",
    shell:
        "poetry run python {input.script}"


rule generate_dashboard:
    input:
        script="scripts/generate_dashboard.py",
        src=SRC_DASHBOARD + SRC_IO + SRC_SUSCEPTIBILITY + SRC_OQMD,
        data="data/processed/synthesis_data.pkl",
        params="results/susceptibility/hc2_fit_parameters.csv",
        chi_data=CHI_DATA_FILES,
        wpf_data=WPF_FILES,
        genai_cifs=GENAI_CIF_FILES,
        xrd_jpegs=XRD_JPG_FILES,
        supercon_cache="data/processed/.supercon_cached",
        ternary_cifs="data/external/oqmd_ternary_phases/.extracted",
        mp_cifs="data/external/mp_ternary_phases/.extracted",
        alex_cifs_pbe="data/external/alexandria_pbe_ternary_phases/.extracted",
        alex_cifs_pbesol="data/external/alexandria_pbesol_ternary_phases/.extracted",
    output:
        index="results/dashboard/index.html",
    shell:
        "poetry run python {input.script}"

# Body-text figure: one representative measured-vs-predicted XRD panel per
# kind of outcome.
rule plot_xrd_representative:
    input:
        script="scripts/plot_xrd_representative.py",
        src=SRC_FIGURES,
        sheets=SAMPLE_SHEETS,
        csv="data/processed/synthesis_data.csv",
        xrd=XRD_FILES,
        cifs=RAW_CIF_FILES,
        oqmd_structures="data/external/oqmd_structures/.extracted",
    output:
        pdf="results/publication_ready/xrd_representative.pdf",
    log:
        "logs/plot_xrd_representative.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


# Provenance of every cached database (versions, dates, checksums): the state the
# paper's numbers were computed from, and the baseline for a later re-pull.
rule record_database_snapshot:
    input:
        script="scripts/record_database_snapshot.py",
        ternary_data=TERNARY_JSON_FILES,
        mp_data=[str(p) for p in Path("data/external/mp_ternary_phases").rglob("*.json")],
        alex_data=[str(p) for p in Path("data/external/alexandria_pbe_ternary_phases").rglob("*.json")]
                 + [str(p) for p in Path("data/external/alexandria_pbesol_ternary_phases").rglob("*.json")],
    output:
        json="data/external/snapshot_provenance.json",
    log:
        "logs/record_database_snapshot.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


# Supplementary table: hull distance in each database and whether its ground
# state is the structure OQMD predicted.
rule make_database_table:
    input:
        script="scripts/make_database_table.py",
        src=["src/synthesizability/paper_samples.py", "src/synthesizability/formula.py"],
        sheets=SAMPLE_SHEETS,
        csv="data/processed/synthesis_data.csv",
        hull="data/processed/hull_distances.csv",
        alex_index="data/external/alexandria_hull_structures/index.csv",
        mp_cifs="data/external/mp_ternary_phases/.extracted",
        oqmd_structures="data/external/oqmd_structures/.extracted",
    output:
        tex="results/publication_ready/database_table.tex",
        csv="results/publication_ready/database_agreement.csv",
    log:
        "logs/make_database_table.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


rule plot_hull_cross_database:
    """OQMD hull distance vs the same quantity in MP and Alexandria, only where
    the database's ground-state structure matches the OQMD prediction."""
    input:
        script="scripts/plot_hull_cross_database.py",
        src="src/synthesizability/formula.py",
        csv="results/publication_ready/database_agreement.csv",
    output:
        pdf="results/publication_ready/hull_cross_database.pdf",
    log:
        "logs/plot_hull_cross_database.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


# Sample table for the paper: formula, hull distance, disorder parameter, outcome.
rule make_sample_table:
    input:
        script="scripts/make_sample_table.py",
        src=["src/synthesizability/paper_samples.py", "src/synthesizability/formula.py"],
        sheets=SAMPLE_SHEETS,
        csv="data/processed/synthesis_data.csv",
        hull="data/processed/hull_distances.csv",
    output:
        tex="results/publication_ready/sample_table.tex",
    log:
        "logs/make_sample_table.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


# Appendix figures: every classified target, one grid per XRD outcome.
rule plot_xrd_comparison_grid:
    input:
        script="scripts/plot_xrd_comparison_grid.py",
        src=SRC_FIGURES,
        sheets=SAMPLE_SHEETS,
        csv="data/processed/synthesis_data.csv",
        xrd=XRD_FILES,
        cifs=RAW_CIF_FILES,
        oqmd_structures="data/external/oqmd_structures/.extracted",
    output:
        pdfs=expand("results/publication_ready/xrd_appendix_{outcome}.pdf", outcome=XRD_OUTCOMES),
    log:
        "logs/plot_xrd_comparison_grid.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"


rule plot_disorder_vs_stability:
    input:
        script="scripts/plot_disorder_vs_stability.py",
        src=["src/synthesizability/paper_samples.py"],
        sheets=SAMPLE_SHEETS,
        csv="data/processed/synthesis_data.csv",
        phonon="data/external/phonon_stability/MANIFEST_w_stability.csv",
    output:
        pdf="results/publication_ready/disorder_vs_stability.pdf",
    log:
        "logs/plot_disorder_vs_stability.log"
    shell:
        "poetry run python {input.script} > {log} 2>&1"
