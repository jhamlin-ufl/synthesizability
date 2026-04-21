# Synthesizability

**Preliminary data repository associated with an upcoming publication.**

This repository contains the analysis pipeline and results for an experimental study of crystal-structure synthesizability: we tested whether materials predicted as stable by high-throughput computational databases actually form into their predicted crystal structures when synthesized by arc melting.

Across numerous intermetallic samples, predicted crystal structures were confirmed by powder X-ray diffraction in only a minority of cases. A disorder-probability model ([Jakob et al. 2025](https://doi.org/10.1002/adma.202514226)) was found to be strongly predictive of synthesis outcome.

A companion paper is in preparation. This repository will be updated to reference the published work when available.

## Exploring the Data

The primary artifact for browsing results is the **interactive dashboard** — a self-contained set of static HTML files that requires no server or installation.

### Download and open the dashboard

```bash
# Download the dashboard (≈45 MB)
curl -L https://github.com/jhamlin-ufl/synthesizability/releases/latest/download/dashboard.zip \
     -o dashboard.zip

unzip dashboard.zip -d synthesizability-dashboard

# Open in your browser
open synthesizability-dashboard/index.html          # macOS
xdg-open synthesizability-dashboard/index.html      # Linux
start synthesizability-dashboard/index.html         # Windows
```

The index page summarizes all samples with sortable/filterable tables. Click any sample row to view its individual page, which includes the XRD pattern, fit, predicted vs. observed phases, thermodynamic stability data, and (where measured) AC susceptibility and superconducting transition data.

## Repository Structure

```
data/
  raw/          raw experimental data (XRD, AC susceptibility)
  processed/    intermediate processed data
  external/     reference databases (OQMD structures, SuperCon, etc.)
results/
  dashboard/    generated static HTML dashboard
  figures/      publication figures
src/synthesizability/   core Python library
scripts/        analysis scripts
Snakefile       pipeline definition (single source of truth)
```

## Reproducing the Analysis

The pipeline is managed by [Snakemake](https://snakemake.readthedocs.io/) and [Poetry](https://python-poetry.org/).

```bash
git clone https://github.com/jhamlin-ufl/synthesizability.git
cd synthesizability
poetry install
poetry run snakemake --cores all
```

Running the full pipeline requires local access to the raw experimental data files (XRD work files and susceptibility data), which are stored in `data/raw/` and are not included in the public repository due to size.

## Citation

A DOI will be added here once the associated paper is published. In the meantime, please cite this repository directly using the Zenodo DOI linked in the GitHub sidebar.

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). You are free to share and adapt the material with appropriate attribution.
