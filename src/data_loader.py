"""Strict CSV validation. User data is never silently coerced or imputed."""
from __future__ import annotations
import csv
import gzip
import io
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
MAX_BYTES = 200 * 1024 * 1024
MAX_CELLS = 10_000_000

def dataframe_fingerprint(df):
    """Hash every row for scientific caching, avoiding sampled-frame hashes."""
    h = hashlib.sha256()
    h.update(json.dumps({"columns": list(map(str, df.columns)), "dtypes": list(map(str, df.dtypes))}, sort_keys=True).encode())
    h.update(pd.util.hash_pandas_object(df, index=True).to_numpy().tobytes())
    return h.hexdigest()

def read_csv(source):
    if isinstance(source, (str, Path)):
        raw = Path(source).read_bytes()
    else:
        raw = source.getvalue() if hasattr(source, "getvalue") else source.read()
    if raw[:2] == b"\x1f\x8b":
        with gzip.GzipFile(fileobj=io.BytesIO(raw)) as f:
            raw = f.read(MAX_BYTES + 1)
    if len(raw) > MAX_BYTES:
        raise ValueError("The expanded CSV exceeds 200 MB. Upload a smaller matrix.")
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError("Save the file as a UTF-8 CSV before uploading.") from exc
    if not text.strip():
        raise ValueError("This CSV is empty.")
    try:
        headers = [v.strip() for v in next(csv.reader(io.StringIO(text)))]
        if any(not h for h in headers) or len(headers) != len(set(headers)):
            raise ValueError("CSV column names must be non-empty and unique.")
        df = pd.read_csv(io.StringIO(text), keep_default_na=False)
    except (pd.errors.ParserError, pd.errors.EmptyDataError, StopIteration) as exc:
        raise ValueError("Could not parse this file. Use a comma-separated CSV with one header row.") from exc
    df.columns = [str(c).strip() for c in df.columns]
    if df.empty:
        raise ValueError("This CSV has a header but no data rows.")
    return df

def validate_expression(df):
    if df.shape[1] < 5:
        raise ValueError("Use Feature/Gene in the first column and at least four sample columns.")
    if (df.shape[1] - 1) * len(df) > MAX_CELLS:
        raise ValueError("This matrix exceeds 10 million values. Select a smaller cohort.")
    ids = df.iloc[:, 0].astype(str).str.strip()
    if ids.eq("").any() or ids.duplicated().any():
        raise ValueError("Feature identifiers must be non-empty and unique. Keep probe IDs when several probes map to one gene.")
    numeric = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError("Expression values must all be finite numbers. Missing values, text, NaN and infinity are not supported.")
    values = numeric.to_numpy(dtype=float)
    if values.min() >= 0 and values.max() > 50 and np.isclose(values, np.round(values)).all():
        raise ValueError("This looks like raw integer counts. GeneDiff requires normalized log2 expression. Analyze raw RNA-seq counts with DESeq2/edgeR first.")
    numeric.index = pd.Index(ids, name="Feature")
    return numeric.astype(float)

def validate_metadata(df, samples):
    if not {"Sample", "Group"}.issubset(df.columns):
        raise ValueError("Metadata needs Sample and Group columns, spelled exactly as shown.")
    df = df.copy()
    for c in ["Sample", "Group"]:
        df[c] = df[c].astype(str).str.strip()
        if df[c].eq("").any():
            raise ValueError(f"Metadata {c} values cannot be empty.")
    if df.Sample.duplicated().any():
        raise ValueError("Each sample must appear only once in the metadata.")
    missing = set(samples) - set(df.Sample)
    extra = set(df.Sample) - set(samples)
    if missing or extra:
        raise ValueError(f"Metadata must match the matrix exactly. Missing samples: {sorted(missing)[:8]}; extra samples: {sorted(extra)[:8]}.")
    return df.set_index("Sample").loc[list(samples)].reset_index()

def validate_amr(df):
    required = {"Isolate_ID", "Bacteria", "Antibiotic", "Result"}
    if not required.issubset(df.columns):
        raise ValueError("AMR CSV needs Isolate_ID, Bacteria, Antibiotic and Result columns.")
    df = df.copy()
    for c in required:
        df[c] = df[c].astype(str).str.strip()
        if df[c].eq("").any():
            raise ValueError(f"{c} cannot contain empty cells.")
    aliases = {"S": "Susceptible", "I": "Intermediate", "R": "Resistant", "SUSCEPTIBLE": "Susceptible", "INTERMEDIATE": "Intermediate", "RESISTANT": "Resistant"}
    df["Result"] = df.Result.str.upper().map(aliases)
    if df.Result.isna().any():
        raise ValueError("Result must be S, I, R, Susceptible, Intermediate or Resistant. Unsupported classifications must be reviewed before upload.")
    if df.duplicated(["Isolate_ID", "Bacteria", "Antibiotic"]).any():
        raise ValueError("Repeated isolate-organism-antibiotic results found. Use globally unique isolate IDs and resolve duplicate tests before upload.")
    if "Year" in df:
        year = pd.to_numeric(df.Year, errors="coerce")
        if year.isna().any() or not np.equal(year, np.floor(year)).all() or not year.between(1900, 2100).all():
            raise ValueError("Year must be a four-digit integer from 1900 to 2100.")
        df["Year"] = year.astype(int)
    else:
        df["Year"] = 0
    for c in ["State", "Host", "Specimen", "Agency"]:
        if c not in df:
            df[c] = "Not reported"
        df[c] = df[c].astype(str).str.strip().replace("", "Not reported")
    # An isolate has one organism and stable source attributes.
    attributes = ["Bacteria", "Year", "State", "Host", "Specimen", "Agency"]
    if df.groupby("Isolate_ID")[attributes].nunique().gt(1).any().any():
        raise ValueError("An Isolate_ID has conflicting organism or source metadata. Make each isolate ID globally unique.")
    return df

def load_gene():
    expression = validate_expression(read_csv(DATA / "gse42568_expression.csv.gz"))
    metadata = validate_metadata(read_csv(DATA / "gse42568_metadata.csv"), expression.columns)
    annotation = read_csv(DATA / "gse42568_annotation.csv.gz")
    return expression, metadata, annotation

def load_amr():
    return validate_amr(read_csv(DATA / "fda_amr.csv.gz"))
