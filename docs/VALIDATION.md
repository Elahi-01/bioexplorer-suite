# Verification performed for this ZIP

The application was exercised with the bundled real datasets and with uploaded
subsets of those real datasets. This report records checks actually run; it
does not claim a live hosted deployment, Docker build or clinical validation.

## Environment

- Python 3.12.14 on Linux.
- Streamlit 1.50.0, pandas 2.2.3, NumPy 2.2.6, SciPy 1.15.3, Plotly 6.3.1,
  openpyxl 3.1.5, pytest 8.4.2.
- Actual browser checks: headless Chromium using Playwright 1.55.0.
- Requirements pin exactly the runtime package versions used for verification.
- Python 3.11 is supplied as the Community Cloud / Docker target; a separate
  local Python 3.11 run was not performed. There is no claim of platform-wide
  compatibility testing on Windows/macOS.

## Automated checks

Command: `python -m pytest -q` from the project directory.

**Result: 19 passed, 0 failed, no reported warnings, in 20.96 seconds.**

Checks cover:

- Welch p-values against direct SciPy calculations and known effect directions.
- Reference/target reversal preserves p-values and reverses effect signs.
- BH values against the documented SciPy example, including missing-p exclusion.
- Constant/no-variance features stay untested rather than becoming false hits.
- FDR cutoff is strict; effect-size cutoff is inclusive.
- Empty/malformed numeric cells, duplicate headers/feature IDs and obvious counts
  are rejected; metadata must match sample IDs.
- Cache keys fingerprint every dataframe row, including changes at row 60,000.
- S/I/R denominators, count-weighted pooling and non-zero interval bounds for
  zero observed resistance.
- Unused categorical combinations never become observed tested groups.
- Duplicate AST keys, unsupported classes and inconsistent source identities
  are rejected.
- Bundled file SHA256 hashes match `manifest.json`; real group sizes and AST
  key uniqueness are validated.
- GeneDiff thresholds, quick cohort and unsupported raw-count UI state.
- AMR filters, multiple organisms, regional view and empty-result state.
- Both standalone app entrypoints and the quick guide.

## Actual browser checks

The Streamlit app was started locally and exercised through Chromium:

| Check | Outcome |
|---|---|
| GeneDiff real study loads and renders charts | Passed |
| GeneDiff heatmap tab opens | Passed |
| AMR real data loads and renders charts | Passed |
| Regional heatmap for one organism | Passed |
| AMR summary download produces an actual CSV | Passed; counts/rates checked |
| Upload a real 8-feature, 3+3 GEO subset + matching metadata | Passed |
| Invalid numeric expression CSV shows validation error | Passed |
| Upload a real 200-row FDA subset | Passed |
| Browser page JavaScript errors during dashboard/export checks | None |
| 390 × 844 mobile viewport horizontal overflow | None |

Screenshots in this folder show the desktop overview, gene heatmap, regional
AMR view and mobile layout. They were visually inspected for readable hierarchy,
contrast, labels and layout. They contain only public dataset information.

## Reproducible default results

These numbers are computed by this app on the bundled snapshot, not copied
from the original study publication.

### GeneDiff

Full GSE42568 cohort, Normal reference, Breast cancer target, FDR < 0.05 and
|log2FC| ≥ 1:

| Metric | Value |
|---|---:|
| Features examined | 54,675 |
| Valid Welch tests / BH family | 54,579 |
| Untested features | 96 |
| Significant features | 6,455 |
| Upregulated in target | 3,409 |
| Downregulated in target | 3,046 |

### AMR Explorer

Default E. coli panel: Ampicillin, Enrofloxacin, Cephalexin, Tetracycline,
Gentamicin and Amoxicillin/Clavulanic Acid; all source years/regions/specimens/
agencies selected:

| Metric | Value |
|---|---:|
| Classified AST results in view | 76,305 |
| Distinct isolate IDs in view | 12,969 |
| Organisms | 1 |
| Antibiotics | 6 |

The entire cleaned FDA bundle contains 499,994 classified results from 27,013
source isolate IDs. These are not unique animal counts.

## Scope of verification

Cloud deployment instructions and files are supplied, but no GitHub repository
or Streamlit Cloud app was created. Docker was not built. Local browser QA and
tests establish that the packaged application runs in the tested environment;
they do not establish biological validity, clinical utility, performance under
large multi-user workloads or the absence of unmodeled study confounders.
