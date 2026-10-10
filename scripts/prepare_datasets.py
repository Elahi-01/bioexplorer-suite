"""Rebuild the bundled real datasets from their primary public sources.

Run from any directory: python scripts/prepare_datasets.py
Use --source-dir to reuse already downloaded sources.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import urllib.request
from datetime import datetime, timezone
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
URLS = {
    "GSE42568_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42568/matrix/GSE42568_series_matrix.txt.gz",
    "GPL570.annot.gz": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    "fda_animal_pathogens.xlsx": "https://www.fda.gov/media/132928/download?attachment",
}

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def build(source_dir):
    source_dir.mkdir(parents=True, exist_ok=True)
    data = ROOT / "data"
    data.mkdir(exist_ok=True)
    for name, url in URLS.items():
        dest = source_dir / name
        if not dest.exists():
            print(f"Downloading {name} ...", flush=True)
            with urllib.request.urlopen(url, timeout=180) as response, open(dest.with_suffix(dest.suffix + '.part'), 'wb') as out:
                while chunk := response.read(1024 * 1024):
                    out.write(chunk)
            dest.with_suffix(dest.suffix + '.part').replace(dest)
    headers = {}
    with gzip.open(source_dir / "GSE42568_series_matrix.txt.gz", "rt") as f:
        for line in f:
            if line.startswith("!Sample_"):
                parts = next(csv.reader([line], delimiter="\t"))
                headers.setdefault(parts[0], parts[1:])
            if line.startswith("!series_matrix_table_begin"):
                expression = pd.read_csv(io.StringIO(f.read().split("!series_matrix_table_end")[0]), sep="\t", index_col=0)
                break
        else:
            raise ValueError("GEO source has no expression table.")
    titles = headers["!Sample_title"]
    samples = headers["!Sample_geo_accession"]
    if list(expression.columns) != samples:
        raise ValueError("GEO sample headers do not match the expression matrix.")
    groups = ["Normal" if t.startswith("Normal breast,") else "Breast cancer" if t.startswith("Breast cancer,") else None for t in titles]
    if None in groups:
        raise ValueError("Unrecognized sample title. Inspect updated source before rebuilding.")
    expression.index.name = "Feature"
    expression.to_csv(data / "gse42568_expression.csv.gz", compression={"method": "gzip", "mtime": 0})
    pd.DataFrame({"Sample": samples, "Group": groups, "Title": titles}).to_csv(data / "gse42568_metadata.csv", index=False)
    with gzip.open(source_dir / "GPL570.annot.gz", "rt") as f:
        content = f.read()
    # GEO annotation format includes a platform table delimited by markers.
    table = content.split("!platform_table_begin\n", 1)[1].split("!platform_table_end", 1)[0]
    annotation = pd.read_csv(io.StringIO(table), sep="\t", dtype=str, keep_default_na=False)
    print("Annotation columns:", list(annotation.columns), flush=True)
    symbol_col = "Gene symbol" if "Gene symbol" in annotation else "Gene Symbol"
    annotation = annotation[["ID", symbol_col]].rename(columns={"ID": "Feature", symbol_col: "Gene_symbol"})
    annotation = annotation[annotation.Feature.isin(expression.index)].drop_duplicates("Feature")
    annotation.to_csv(data / "gse42568_annotation.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    print(f"GEO: {expression.shape[0]:,} probe sets, {expression.shape[1]} samples.", flush=True)
    print("Reading the FDA workbook ...", flush=True)
    raw = pd.read_excel(source_dir / "fda_animal_pathogens.xlsx", sheet_name=0)
    print("FDA rows:", len(raw), "interpretations:", raw.Interpretation.value_counts(dropna=False).to_dict(), flush=True)
    print("Reported organisms:", raw.Genus.value_counts(dropna=False).to_dict(), flush=True)
    rename = {"Agency": "Agency", "Sample ID": "Original_sample_ID", "State": "State", "Host Species": "Host", "Collection Source": "Specimen", "Genus": "Bacteria", "Drug Name": "Antibiotic", "MIC": "MIC", "Sign": "MIC_sign", "Interpretation": "Result", "Year": "Year"}
    amr = raw[list(rename)].rename(columns=rename)
    amr["Result"] = amr.Result.astype("string").str.strip()
    allowed = ["Susceptible", "Intermediate", "Resistant"]
    known = amr.Result.isin(allowed)
    excluded = amr.loc[~known, "Result"].fillna("Missing").value_counts().to_dict()
    amr = amr.loc[known].copy()
    required = ["Agency", "Original_sample_ID", "Bacteria", "Antibiotic", "Year", "Result"]
    missing = int(amr[required].isna().any(axis=1).sum())
    amr = amr.dropna(subset=required)
    amr["Year"] = amr.Year.astype(int)
    amr["Bacteria"] = amr.Bacteria.replace({"E. coli": "Escherichia coli", "S. pseudintermedius": "Staphylococcus pseudintermedius"})
    for col in ["State", "Host", "Specimen"]:
        amr[col] = amr[col].fillna("Not reported").astype(str).str.strip()
    amr["Isolate_ID"] = amr.Agency.astype(str) + "|" + amr.Year.astype(str) + "|" + amr.Original_sample_ID.astype(str)
    attributes = ["Bacteria", "Year", "State", "Host", "Specimen", "Agency"]
    meta_counts = amr.groupby("Isolate_ID")[attributes].nunique()
    ambiguous_ids = meta_counts.index[meta_counts.gt(1).any(axis=1)]
    ambiguous = amr.Isolate_ID.isin(ambiguous_ids)
    ambiguous_rows = int(ambiguous.sum())
    if ambiguous_rows:
        amr.loc[ambiguous].to_csv(data / "fda_ambiguous_metadata_excluded.csv", index=False)
        amr = amr.loc[~ambiguous].copy()
    else:
        (data / "fda_ambiguous_metadata_excluded.csv").unlink(missing_ok=True)
    key = ["Isolate_ID", "Bacteria", "Antibiotic"]
    # Do not silently deduplicate repeated observations with different interpretations.
    repeated = amr.duplicated(key, keep=False)
    conflicting = amr.loc[repeated].groupby(key).Result.nunique()
    conflict_keys = conflicting[conflicting.gt(1)].reset_index()[key]
    if len(conflict_keys):
        audit = amr.merge(conflict_keys, on=key)
        audit.to_csv(data / "fda_conflicting_results_excluded.csv", index=False)
        marked = amr.merge(conflict_keys.assign(_conflict=True), on=key, how="left")
        conflict_rows = int(marked._conflict.notna().sum())
        amr = marked.loc[marked._conflict.isna()].drop(columns="_conflict")
    else:
        conflict_rows = 0
        (data / "fda_conflicting_results_excluded.csv").unlink(missing_ok=True)
    duplicates = int(amr.duplicated(key).sum())
    amr = amr.drop_duplicates(key)
    cols = ["Isolate_ID", "Bacteria", "Antibiotic", "Result", "Year", "State", "Host", "Specimen", "Agency", "MIC", "MIC_sign", "Original_sample_ID"]
    amr[cols].to_csv(data / "fda_amr.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    files = sorted(data.glob("*.csv*"))
    manifest = {
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {name: {"url": url, "sha256": sha256(source_dir / name)} for name, url in URLS.items()},
        "gene": {"accession": "GSE42568", "platform": "GPL570", "features": len(expression), "samples": len(samples), "groups": pd.Series(groups).value_counts().to_dict(), "scale": "Submitter-provided log2 GC-RMA signal intensity", "selection": "All probe sets and all samples; no rounding, filtering, re-normalization or re-logging."},
        "amr": {"source_rows": len(raw), "bundled_ast_rows": len(amr), "isolates": amr.Isolate_ID.nunique(), "years": sorted(amr.Year.unique().tolist()), "organisms": sorted(amr.Bacteria.unique().tolist()), "regions": sorted(amr.State.unique().tolist()), "interpretations_excluded": excluded, "missing_key_rows_excluded": missing, "ambiguous_metadata_isolates_excluded": len(ambiguous_ids), "ambiguous_metadata_rows_excluded": ambiguous_rows, "conflicting_key_count": len(conflict_keys), "conflicting_result_rows_excluded": conflict_rows, "same_key_same_result_duplicates_removed": duplicates, "result_counts": amr.Result.value_counts().to_dict(), "selection": "All workbook years and organisms; retain source-reported S/I/R only. Exclude isolates with conflicting source metadata; exclude all rows for conflicting AST keys; collapse same-result duplicate keys. No MIC breakpoint reinterpretation."},
        "files": {p.name: {"sha256": sha256(p), "bytes": p.stat().st_size} for p in files},
    }
    (data / "manifest.json").write_text(json.dumps(manifest, indent=2, default=int) + "\n")
    print(json.dumps(manifest["amr"], indent=2, default=int), flush=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=ROOT / ".source_cache")
    build(parser.parse_args().source_dir)
