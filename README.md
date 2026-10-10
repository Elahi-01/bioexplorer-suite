# BioExplorer: GeneDiff + AMR Explorer

Two deployable Streamlit research apps in one repository. The default workspace
opens GeneDiff; switch to AMR Explorer or the quick guide in the sidebar.

Both apps ship with **real, source-linked public data**. There is no generated
demo data disguised as real data, no login requirement, and no API key. After
installation, the bundled analyses work without downloading anything at runtime.


## What is included

| App | Bundled source | Main features |
|---|---|---|
| GeneDiff | NCBI GEO GSE42568; normalized log2 GC-RMA microarray data | Independent Welch tests, BH FDR, volcano plot, clustered feature heatmap, feature search/box plot, sample PCA, expression QC, CSV/JSON/HTML exports |
| AMR Explorer | FDA Animal Pathogen AMR workbook; canine Vet-LIRN/NAHLN surveillance, 2017–2024 | Organism/drug/year/region/specimen/agency filters, resistant/tested counts, Wilson intervals, annual trends, S/I/R mix, regional heatmaps, CSV/JSON/HTML exports |

GeneDiff is an expression research explorer, not a disease diagnosis system.
AMR Explorer is descriptive surveillance software, not a prescribing system.

## Run locally

Use Python **3.11 or 3.12**. Tested dependency versions are pinned in
`requirements.txt`. Commands below assume the current directory contains
`app.py`, `requirements.txt` and `data/`.

### Windows

Double-click `run_windows.bat`, or run:

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe -m streamlit run app.py
```

No PowerShell activation-policy change is needed because these commands call
the virtual environment's Python directly.

### macOS / Linux

```bash
bash run_mac_linux.sh
```

Or:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m streamlit run app.py
```

Open the local URL printed in the terminal, normally `http://localhost:8501`.
The initial run loads the bundled data; later runs reuse cached data/results.

## Publish with GitHub + Streamlit Community Cloud

1. Extract the ZIP and create a GitHub repository.
2. Upload **the contents of `bioexplorer_suite/`** to the repository root,
   including `src/`, `data/`, `.streamlit/`, `app.py` and `requirements.txt`.
   Do not upload `.venv/`, source-download caches or the ZIP as the only file.
3. Sign in at [Streamlit Community Cloud](https://share.streamlit.io/) with your
   GitHub account. Choose **Create app** and the repository/branch.
4. Set the main file path to **`app.py`**. In Advanced settings choose
   **Python 3.11**. No secrets are needed.
5. Deploy. The generated `*.streamlit.app` URL is the live app you can share.

The app files and bundled data are below GitHub's ordinary 100 MB per-file
limit; Git LFS is not required. Keep the compressed `.csv.gz` data compressed.

**Two separate URLs:** deploy the same repository twice, choosing
`genediff_app.py` for GeneDiff and `amr_app.py` for AMR Explorer.

**Nested repository folder:** if you retain `bioexplorer_suite/` in the GitHub
repository, use `bioexplorer_suite/app.py` as the main file path and retain
`requirements.txt` alongside that entrypoint. Putting this project's contents
at the repository root is the simplest option.

GitHub stores the source code. **GitHub Pages does not run Python/Streamlit**;
use Streamlit Community Cloud or another Python host for a live analysis app.
No hosted deployment has been performed as part of this ZIP deliverable.

Official instructions:
[Streamlit deployment guide](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app).

## Upload formats

### GeneDiff: expression CSV + optional metadata CSV

The first column contains a unique feature ID. Every other column is a sample
with finite numeric **normalized log2 expression**. CSV and gzip CSV are accepted.

```csv
Feature,H1,H2,H3,D1,D2,D3
TP53,5.1,5.3,5.0,7.2,7.5,7.1
EGFR,3.2,3.4,3.3,6.1,6.3,6.0
BRCA1,7.1,7.0,7.3,4.4,4.5,4.3
ACTB,8.2,8.1,8.3,8.2,8.0,8.1
```

These README numbers are illustrative only. The bundled GEO data is real.

```csv
Sample,Group
H1,Healthy
H2,Healthy
H3,Healthy
D1,Disease
D2,Disease
D3,Disease
```

Metadata must match the expression sample names exactly. Without metadata,
select reference and target samples manually. Unassigned samples are excluded
and counted explicitly. With more than two groups, choose a two-group comparison.
At least two independent samples per group are required; three or more are
preferable. Confirm independent biological samples for uploaded matrices.

Missing values, duplicate feature/sample identifiers, text expression cells,
infinities and obvious raw-count matrices are rejected. No silent imputation,
automatic re-normalization or automatic logarithm is applied. Scale cannot be
reliably inferred from numbers alone; users must identify it correctly.

**Raw RNA-seq counts are not supported.** Analyze them using an appropriate
DESeq2/edgeR pipeline, including normalization, a suitable design and count-aware
inference. Simply taking log2 of raw counts does not make this Welch workflow
appropriate. Paired samples or important confounding variables need another model.

### AMR Explorer: one row per classified AST result

Required column names are exact:

```csv
Isolate_ID,Bacteria,Antibiotic,Result,Year,State,Host,Specimen,Agency
ISO-001,Escherichia coli,Ciprofloxacin,R,2024,Region A,Human,Blood,Study A
ISO-002,Escherichia coli,Ciprofloxacin,S,2024,Region A,Human,Blood,Study A
ISO-003,Escherichia coli,Ciprofloxacin,I,2024,Region B,Human,Blood,Study A
```

These format examples are illustrative. Required: `Isolate_ID`, `Bacteria`,
`Antibiotic`, `Result`. Optional: `Year`, `State`, `Host`, `Specimen`, `Agency`,
and additional informational fields such as `MIC`.

`Result` accepts S/I/R or Susceptible/Intermediate/Resistant, case-insensitively.
The I code remains separate and retains its source-defined meaning; under EUCAST,
I means susceptible with increased exposure. The app does not apply breakpoints
to MIC values or merge I into resistance. For source-specific alternate labels,
prepare an explicitly documented mapping before upload.

Each isolate ID must be globally unique and must have consistent organism,
year, region and source metadata. Duplicate isolate–organism–antibiotic rows are
rejected, including identical duplicates, so users explicitly resolve them.
Do not upload personal identifiers or private patient information to a public app.

## Interpret the analyses

### GeneDiff

- `log2FC = target mean log2 expression − reference mean log2 expression`.
  This estimates the log2 ratio of geometric means. +1 is approximately 2×;
  −1 is approximately ½×. It is not a ratio of arithmetic means of log values.
- P-values come from a two-sided, independent-group Welch t-test.
- A feature whose two groups both have variance ≤ 1e−12 is not testable. Its
  effect remains available, while p/FDR stay missing and it cannot be significant.
- BH correction uses every finite test p-value in the selected comparison,
  before filtering by effect size or search text.
- Default significance is **FDR < 0.05 AND |log2FC| ≥ 1**. Cutoffs are adjustable.
- FDR is a multiple-testing criterion, not a posterior probability for an
  individual feature, proof of causation or a universal biological cutoff.
- Heatmap rows are significant features ranked by FDR and then hierarchically
  ordered using average-linkage Euclidean distances on row z-scores. If none pass,
  the most variable features are displayed and labeled as exploratory.
- Heatmap colors use each feature's z-score across the compared samples, clipped
  to ±3 for display. Testing still uses the original normalized log2 values.
- PCA centers up to 2,000 most variable features and uses NumPy SVD. It does not
  use labels to select features, batch-correct, or classify future samples.
- Microarray probes mapping to several symbols retain their source annotation.
  Probes are not collapsed into unique genes; feature counts are not gene counts.

This prototype has no covariate adjustment, empirical Bayes variance moderation,
paired testing, pathway enrichment or count-aware RNA-seq analysis. Consider
design-aware [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
for publication-oriented microarray analysis. The independent Welch approach is
an exploratory implementation of the requested prototype.

### AMR Explorer

- Each organism–antibiotic percentage is `100 × resistant / (S + I + R)`.
- Untested and non-interpretable observations are not assumed susceptible.
- Groups are pooled using counts, not averages of subgroup percentages.
- There is no overall resistance KPI averaged across different antibiotics.
- One isolate can contribute to multiple antibiotic denominators. Isolate IDs
  are not unique-animal or unique-patient IDs.
- Wilson 95% intervals describe binomial sample uncertainty. They do not adjust
  for surveillance selection, laboratory clustering, repeat animals, or changes
  in testing/breakpoints. They are not population prevalence intervals.
- Minimum tested count defaults to 30 for charts. Small groups are suppressed
  visually but retained in exports. Missing/suppressed heatmap cells are blank;
  absent years are gaps in trend lines.
- Geographic comparisons hold the organism fixed. Source regions include
  U.S. states and Canada; the snapshot is not a representative human survey.
- The FDA workbook reports S and R only. A zero intermediate bar means no
  Intermediate-classified records were supplied, not biological absence.

Full data sources and cleaning details: [DATA_SOURCES.md](DATA_SOURCES.md).

## Exports

GeneDiff exports full, significant and filtered result CSVs, analysis settings
JSON, an interactive volcano HTML, an interactive heatmap HTML and input data.
AMR Explorer exports all pair summaries, filtered AST rows, filter settings JSON,
an interactive resistance HTML and the full bundled compressed dataset.

HTML plot exports embed Plotly and work offline. Use the camera icon in a plot's
toolbar to save a PNG. CSV exports neutralize formula-like text cells for safer
spreadsheet opening; numeric results remain numeric.

## Project files

| Path | Purpose |
|---|---|
| `app.py` | Combined workspace / default deployment entrypoint |
| `genediff_app.py`, `amr_app.py` | Optional standalone deployments |
| `src/data_loader.py` | CSV and scientific input validation |
| `src/differential.py` | Welch analysis, FDR, PCA |
| `src/amr_analysis.py` | S/I/R counts and Wilson intervals |
| `src/plots.py`, `src/ui.py` | Plotly figures and responsive styling |
| `src/gene_page.py`, `src/amr_page.py`, `src/guide_page.py` | Streamlit interfaces |
| `data/` | Real prepared CSVs, annotation, audit CSVs and SHA256 manifest |
| `scripts/prepare_datasets.py` | Re-download and rebuild from primary sources |
| `tests/` | Numerical, validation, provenance and app interaction tests |
| `.streamlit/config.toml` | Theme, upload cap and hosting configuration |
| `.github/workflows/test.yml` | GitHub Actions checks |
| `Dockerfile` | Optional container deployment |
| `docs/` | Verified interface screenshots and validation notes |

## Run checks or refresh datasets

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
python scripts/prepare_datasets.py
```

Refreshing downloads the original GEO matrix, platform annotation and FDA
workbook into `.source_cache/`. The script rewrites the prepared datasets and
manifest, including source and output SHA256 values. Source files can change;
review new exclusions and group counts before committing a refreshed snapshot.

Optional Docker:

```bash
docker build -t bioexplorer .
docker run --rm -p 8501:8501 bioexplorer
```

The Docker recipe is supplied but a Docker build was not run during artifact
verification. See `docs/VALIDATION.md` for the checks actually performed.

## Troubleshooting

- **Module not found:** install `requirements.txt` using the same Python that
  launches Streamlit; use `python -m streamlit`, not a different global executable.
- **No data on Community Cloud:** upload the `data/` folder and retain its exact
  filenames. No external data API is called by the app at runtime.
- **Large initial load:** let the first run finish, or use GeneDiff's 10+10 cohort
  and narrower AMR filters. Cached results are bounded to keep memory manageable.
- **No significant genes:** inspect QC, sample design and effect sizes. Do not
  treat changing the thresholds until a feature appears as confirmatory evidence.
- **No AMR rows:** select at least one organism/drug and widen filters. Some
  organism–drug combinations were never tested/classified in the source.
- **Only one AMR year:** trends require at least two reported years.

## License and attribution

Generated application code is provided under the MIT license in `LICENSE`.
Bundled third-party datasets and annotations retain their source terms; the
code license does not relicense them. Keep `DATA_SOURCES.md`, `data/manifest.json`
and the source citations when redistributing the project or presenting results.
