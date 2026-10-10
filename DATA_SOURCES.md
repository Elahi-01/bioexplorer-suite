# Dataset provenance and transformations

This project bundles prepared **real public observations**, not simulated demo
data. Primary download URLs, acquisition time, original-file SHA256 hashes and
prepared-file SHA256 hashes are in `data/manifest.json`. The original large
source downloads are not included in the ZIP; `scripts/prepare_datasets.py`
downloads and rebuilds them.

## GeneDiff: GSE42568

| Item | Details |
|---|---|
| Repository | NCBI Gene Expression Omnibus (GEO) |
| Study accession | GSE42568, Breast Cancer Gene Expression Analysis |
| Platform | GPL570, Affymetrix Human Genome U133 Plus 2.0 |
| Assay | Expression profiling by array, not RNA-seq |
| Expression scale | Submitter's log2 GC-RMA signal intensity |
| Samples | 17 normal breast + 104 breast-cancer biopsies = 121 |
| Rows | 54,675 probe sets, not 54,675 unique genes |
| Default cohort | All source samples and features |
| Quick cohort | First 10 Normal and first 10 Breast cancer samples in original source order; all features retained |
| Prepared files | `gse42568_expression.csv.gz`, `gse42568_metadata.csv`, `gse42568_annotation.csv.gz` |

Primary sources:

- [Study accession and study-linked publication](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568).
- [Original series matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42568/matrix/GSE42568_series_matrix.txt.gz).
- [Example sample explicitly describing the log2 GC-RMA values](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1045202).
- [Platform annotation](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz).

The preparation script extracts the GEO matrix without rounding, applying a
logarithm, re-normalizing, selecting significant rows or discarding samples.
The expression row ID column is renamed `Feature`; sample accession IDs remain
unchanged. `Group` is assigned from the original sample title prefix:
`Normal breast,` → `Normal`; `Breast cancer,` → `Breast cancer`. Original sample
titles are retained in metadata. The build stops on an unrecognized title.

Probe-to-symbol annotations use GEO's GPL570 annotation table. The script keeps
`ID` → `Feature` and `Gene symbol` → `Gene_symbol`, restricted to matrix probe IDs.
Multi-symbol mappings stay multi-symbol; no transcript/gene aggregation is done.
Source annotation age and ambiguity should be considered before interpreting a
particular probe. The app keeps the probe ID visible even when a symbol is shown.

Expression observations remain unchanged. QC and plotting create temporary
derived representations; they do not overwrite the bundled matrix. Welch
results are new exploratory analyses and need not reproduce the original paper.

Suggested dataset attribution: NCBI GEO accession **GSE42568**, submitted by
Clarke C and colleagues; cite the study-linked publication from the accession
page for academic work. Also identify GPL570 when citing probe annotation.

## AMR Explorer: FDA Animal Pathogen AMR data

| Item | Details |
|---|---|
| Publisher | U.S. Food and Drug Administration |
| Program | Animal Pathogen AMR data, Vet-LIRN and USDA NAHLN |
| Original workbook | FDA media 132928, sheet `2017-2024_data` |
| Population | Clinical bacterial isolates from dogs |
| Source location labels | U.S. states and Canada, plus Not reported |
| Reported organisms | E. coli; S. pseudintermedius; Klebsiella spp; Salmonella |
| Source rows | 533,831 AST-result rows |
| Prepared classified rows | 499,994 |
| Prepared distinct isolate IDs | 27,013 |
| Years | 2017–2024 |
| Source reported classes | Susceptible, Resistant, Non-Interpretable |
| Prepared classes | 363,535 Susceptible; 136,459 Resistant; no Intermediate rows supplied |
| Prepared file | `fda_amr.csv.gz` |

Primary sources:

- [Original FDA workbook](https://www.fda.gov/media/132928/download?attachment).
- [FDA source index and attribution guidance](https://www.fda.gov/animal-veterinary/national-antimicrobial-resistance-monitoring-system/integrated-reportssummaries).
- [Vet-LIRN program](https://www.fda.gov/animal-veterinary/science-research/veterinary-laboratory-investigation-and-response-network).

This is a historical source snapshot, not a claim that it is the most recent
available surveillance year. The app is not a live API client. The source may
be updated; the bundled hashes identify the exact acquisition used here.

### Preparation rules, in order

1. Select source columns needed for identity, organism, drug, interpretation,
   year, region, host, specimen, reporting agency and MIC display.
2. Retain source-reported S/I/R classes. Exclude **33,412 Non-Interpretable** rows.
   Nothing is reclassified from MIC values. The current workbook has no I rows.
3. Check required key fields. **0 rows** are removed for missing keys.
4. Preserve organism resolution: E. coli → Escherichia coli and
   S. pseudintermedius → Staphylococcus pseudintermedius. Klebsiella spp is not
   relabeled as K. pneumoniae; Salmonella is not assigned an unreported species.
5. Fill empty region/host/specimen labels with `Not reported`. Construct a global
   app isolate ID as `Agency|Year|Original_sample_ID`. The original source ID is
   retained. This identifies a source isolate, not a unique animal.
6. Exclude **422 classified rows from 23 isolate IDs** whose organism/source
   metadata conflicts across antibiotic records. No location/specimen is guessed
   or repaired. Audit: `fda_ambiguous_metadata_excluded.csv`.
7. Check the `Isolate_ID + Bacteria + Antibiotic` key. Exclude **both rows** for
   **1 key** with conflicting Susceptible/Resistant results. Audit:
   `fda_conflicting_results_excluded.csv`. The other drug results from that
   isolate remain usable.
8. Remove **1 same-key, same-result duplicate**. Then validate uniqueness and
   consistent isolate metadata in the app's loader.

Accounting: `533,831 − 33,412 − 422 − 2 − 1 = 499,994` prepared rows.
Row-level audit CSVs are public source observations and contain no new invented
results. The full non-interpretable rows can be recovered from the original
workbook using the preparation script/source URL.

Source field mappings:

| Original | Prepared |
|---|---|
| Sample ID | Original_sample_ID, plus constructed Isolate_ID |
| Genus | Bacteria |
| Drug Name | Antibiotic |
| Interpretation | Result |
| Host Species | Host |
| Collection Source | Specimen |
| Sign | MIC_sign |
| Agency, State, Year, MIC | Retained |

FDA's `Genus` column contains labels at different taxonomic resolutions. The app
uses these reported labels and does not infer species. Classifications depend
on the original surveillance testing rules; retaining a source label does not
establish current CLSI/EUCAST clinical appropriateness for a different host or
specimen. FDA source labels S/R should not be assumed to provide a complete
clinical S/I/R distinction.

Suggested attribution: **Food and Drug Administration, Animal Pathogen AMR
Data (Vet-LIRN/NAHLN), via the NARMS data download index**, 2017–2024 snapshot;
include the original URL, acquisition date from `manifest.json`, and these
preparation rules when reporting app-derived estimates.

## Statistical and deployment references

- [SciPy independent t-test / Welch option](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html).
- [SciPy Benjamini–Hochberg FDR adjustment](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html).
- [Bioconductor limma](https://bioconductor.org/packages/release/bioc/html/limma.html).
- [Bioconductor DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
- [EUCAST S/I/R definitions](https://www.eucast.org/bacteria/clinical-breakpoints-and-interpretation/definition-of-s-i-and-r/).
- [Official Streamlit Community Cloud deployment](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app).

The original screenshots were layout references only. The apps are independently
implemented and do not copy QIAGEN branding, connect to IPA, or claim affiliation
with QIAGEN, NCBI, FDA or the reference dashboards.

## Terms and redistribution

The application MIT license covers generated code. The source datasets and
annotations remain subject to their original source terms. Public availability
is documented here; the project does not invent a Creative Commons license or
apply its code license to third-party data. Retain provenance and source credits.
