"""
SV Clonal Architecture Analysis
================================
Classifies structural variants as:
  - truncal  : present in ALL tumor samples of a patient
  - shared   : present in 2+ but not all tumor samples
  - private  : present in exactly 1 tumor sample

Also computes:
  - Tumor VAF  = VCF_TALT / (VCF_TALT + VCF_TREF)
  - Normal VAF = VCF_NALT / (VCF_NALT + VCF_NREF)

Outputs:
  1. TSV with classification + VAFs appended

Usage:
  python sv_clonal_analysis.py --input svs.tsv --output results/ --bp-tolerance 3

Parameters:
  --bp-tolerance   Max breakpoint distance (bp) to consider two SVs the same event.
                   Set to 0 to require exact coordinate match. We accepted a missmatch of 3 bp. Default: 50.
"""

import argparse
import os
import sys
import warnings
import numpy as np
import pandas as pd


warnings.filterwarnings("ignore")


# CLI

def parse_args():
    p = argparse.ArgumentParser(description="SV clonal architecture analysis")
    p.add_argument("--input",  "-i", required=True,  help="Input SV TSV file")
    p.add_argument("--output", "-o", default="sv_results", help="Output directory")
    p.add_argument("--bp-tolerance", "-t", type=int, default=50,
                   help="Max bp distance to consider two SVs identical (default: 50). "
                        "Use 0 for exact match.")
    p.add_argument("--sep", default="\t", help="Column separator (default: tab)")
    return p.parse_args()



# LOAD & VALIDATE

REQUIRED_COLS = ["chr1", "pos1", "chr2", "pos2", "str1", "str2", "class",
                 "Normal_ID", "Tumor_ID",
                 "VCF_TALT", "VCF_TREF", "VCF_NALT", "VCF_NREF"]

def load_data(path: str, sep: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep, dtype=str, low_memory=False)
    # Normalise column names (strip whitespace, lowercase)
    df.columns = df.columns.str.strip()

    # Tolerate mixed case
    col_map = {c.lower(): c for c in df.columns}
    rename = {}
    for req in REQUIRED_COLS:
        if req not in df.columns and req.lower() in col_map:
            rename[col_map[req.lower()]] = req
    if rename:
        df = df.rename(columns=rename)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: Missing required columns: {missing}\n"
                 f"Found columns: {list(df.columns)}")

    # Coerce numeric columns
    for col in ["pos1", "pos2", "VCF_TALT", "VCF_TREF", "VCF_NALT", "VCF_NREF"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    print(f"Loaded {len(df)} SVs across "
          f"{df['Normal_ID'].nunique()} patients / "
          f"{df['Tumor_ID'].nunique()} tumor samples")
    return df


# VAF CALCULATION

def compute_vaf(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    # Tumor VAF
    t_total = df["VCF_TALT"] + df["VCF_TREF"]
    df["Tumor_VAF"] = np.where(t_total > 0, df["VCF_TALT"] / t_total, np.nan)
    # Normal VAF
    n_total = df["VCF_NALT"] + df["VCF_NREF"]
    df["Normal_VAF"] = np.where(n_total > 0, df["VCF_NALT"] / n_total, np.nan)
    return df


# SV MATCHING

def svs_match(sv1: pd.Series, sv2: pd.Series, tol: int) -> bool:
    """
    Two SVs are considered the same event if:
      - Same chromosomes (chr1, chr2)
      - Same strands (str1, str2)
      - Same SV class
      - Both breakpoints within `tol` bp of each other
    When tol=0 coordinates must be identical.
    """
    if sv1["class"] != sv2["class"]:
        return False
    if sv1["chr1"] != sv2["chr1"] or sv1["chr2"] != sv2["chr2"]:
        return False
    if sv1["str1"] != sv2["str1"] or sv1["str2"] != sv2["str2"]:
        return False
    # For inter-chromosomal SVs the strands already disambiguate orientation
    d1 = abs(sv1["pos1"] - sv2["pos1"])
    d2 = abs(sv1["pos2"] - sv2["pos2"])
    return d1 <= tol and d2 <= tol


def cluster_svs_within_patient(patient_df: pd.DataFrame, tol: int) -> pd.Series:
    """
    Greedy single-linkage clustering of SVs within one patient.
    Returns a Series (same index as patient_df) mapping each row → cluster_id.
    """
    n = len(patient_df)
    rows = patient_df.reset_index(drop=False)  # keeps original index
    parent = list(range(n))  # union-find

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for i in range(n):
        for j in range(i + 1, n):
            if svs_match(rows.iloc[i], rows.iloc[j], tol):
                union(i, j)

    cluster_ids = [find(i) for i in range(n)]
    return pd.Series(cluster_ids, index=rows["index"])


# CLASSIFICATION

def classify_svs(df: pd.DataFrame, tol: int) -> pd.DataFrame:
    """
    For each patient (Normal_ID):
      1. Cluster SVs across all tumor samples
      2. Count how many distinct Tumor_IDs share each cluster
      3. Assign truncal / shared / private
    """
    df = df.copy()
    df["cluster_id"] = pd.NA
    df["n_samples_sharing"] = pd.NA
    df["clonality"] = pd.NA
    df["shared_by_tumors"] = pd.NA

    for patient_id, patient_df in df.groupby("Normal_ID"):
        tumor_ids = patient_df["Tumor_ID"].unique()
        n_tumors = len(tumor_ids)

        if n_tumors == 1:
            # Only one sample — everything is "private" by definition
            df.loc[patient_df.index, "cluster_id"] = range(len(patient_df))
            df.loc[patient_df.index, "n_samples_sharing"] = 1
            df.loc[patient_df.index, "clonality"] = "private"
            df.loc[patient_df.index, "shared_by_tumors"] = tumor_ids[0]
            continue

        # Cluster
        cluster_series = cluster_svs_within_patient(patient_df, tol)
        df.loc[patient_df.index, "cluster_id"] = cluster_series.values

        # Annotate globally unique cluster IDs (patient prefix)
        df.loc[patient_df.index, "cluster_id"] = (
            patient_id + "_" + cluster_series.astype(str)
        )

        # For each cluster: which tumor samples are represented?
        temp = patient_df.copy()
        temp["_cluster"] = cluster_series.values
        cluster_tumor_map = (
            temp.groupby("_cluster")["Tumor_ID"]
            .apply(lambda x: sorted(set(x)))
            .to_dict()
        )

        for idx, row in patient_df.iterrows():
            # find by matching the cluster series value
            matching_clusters = [k for k, v in cluster_tumor_map.items()
                                  if k == cluster_series[idx]]
            if matching_clusters:
                tumors_sharing = cluster_tumor_map[matching_clusters[0]]
            else:
                tumors_sharing = [row["Tumor_ID"]]

            n_sharing = len(tumors_sharing)
            df.at[idx, "n_samples_sharing"] = n_sharing
            df.at[idx, "shared_by_tumors"] = ";".join(tumors_sharing)

            if n_sharing == n_tumors:
                df.at[idx, "clonality"] = "truncal"
            elif n_sharing >= 2:
                df.at[idx, "clonality"] = "shared"
            else:
                df.at[idx, "clonality"] = "private"

    return df

# MAIN

def main():
    args = parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Load
    df = load_data(args.input, args.sep)

    # VAF
    print("Computing VAFs …")
    df = compute_vaf(df)

    # Classify
    print(f"Classifying SVs (breakpoint tolerance = {args.bp_tolerance} bp) …")
    df = classify_svs(df, tol=args.bp_tolerance)

    # Write annotated TSV
    out_tsv = os.path.join(args.output, "svs_classified.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)

    print(f"\nAnnotated SV table → {out_tsv}")
    print("Done.")


if __name__ == "__main__":
    main()