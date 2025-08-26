#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import argparse
import json
from pathlib import Path

def load_patient_json(patient_json_path):
    with open(patient_json_path) as f:
        return json.load(f)

def read_maf(maf_path):
    return pd.read_csv(maf_path, sep='\t', low_memory=False)

def select_maf_file_and_allele_cols(sample, genotyped_mafs):
    """
    Identify the genotyped MAF file corresponding to this sample.
    """
    sample_id = sample["sample_id"]
    assay_type = sample["assay_type"]
    tumor_normal = sample["tumor_normal"]

    # Determine suffix and bam_type
    if assay_type in ["research_access", "clinical_access"]:
        if tumor_normal == "tumor":
            suffix = "-SIMPLEX-DUPLEX_genotyped.maf"
            bam_type = "simplex-duplex"
            alt_count_col = "t_alt_count_fragment_simplex_duplex"
            total_count_col = "t_total_count_fragment_simplex_duplex"
        else:
            suffix = "-STANDARD_genotyped.maf"
            bam_type = "unfiltered"
            alt_count_col = "t_alt_count_fragment"
            total_count_col = "t_total_count_fragment"
    elif assay_type == "clinical_impact":
        suffix = "-STANDARD_genotyped.maf"
        bam_type = "standard"
        alt_count_col = "t_alt_count_fragment"
        total_count_col = "t_total_count_fragment"
    else:
        return None, None, None, None

    # Look for matching file in provided list
    for maf_file in genotyped_mafs:
        if maf_file.endswith(f"{sample_id}{suffix}"):
            return maf_file, bam_type, alt_count_col, total_count_col

    return None, None, None, None

def get_variant_counts(maf_df, alt_count_col, total_count_col):
    df = maf_df.copy()
    df["alt_count"] = df[alt_count_col]
    df["total_count"] = df[total_count_col]
    df["vaf"] = df["alt_count"] / df["total_count"].replace(0, np.nan)
    return df

    # Add call_status column

def get_call_status(row, access_min_cov=100, impact_min_cov=50):
    """Determine call status based on assay-specific coverage thresholds."""
    sample_id = str(row["sample_id"]) if not pd.isna(row["sample_id"]) else ""
    called_in = str(row["Called_In"]) if not pd.isna(row["Called_In"]) else ""
    vaf = row.get("vaf", np.nan)
    total_count = row.get("total_count", np.nan)
    assay = row.get("assay", "")
    
    # Select threshold based on assay type
    min_cov = impact_min_cov if assay == "IMPACT" else access_min_cov
    
    if sample_id and sample_id in [x.strip() for x in called_in.split(";")]:
        return "called"
    elif not pd.isna(total_count) and total_count < min_cov:
        return "low_coverage"
    elif not pd.isna(vaf) and vaf > 0:
        return "genotyped"
    else:
        return ""


def aggregate_variants(patient_json_path, genotyped_mafs, union_calls_maf, output_path, access_min_cov=100, impact_min_cov=50):
    # Load patient JSON
    patient_data = load_patient_json(patient_json_path)
    patient_id = patient_data["combined_id"]
    union_df = pd.read_csv(union_calls_maf, sep="\t", low_memory=False)
    union_df_cols = [
        "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", 
        "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", 
        "HGVSp", "HGVSp_Short", "HGVSc", "Clinical", "Called_In"]
    union_df = union_df[union_df_cols]
    # Ensure consistent types for merge
    union_df["Chromosome"] = union_df["Chromosome"].astype(str)
    union_df["Start_Position"] = union_df["Start_Position"].astype(int)
    union_df["End_Position"] = union_df["End_Position"].astype(int)

    combined_list = []

    for sample_name, sample in patient_data["samples"].items():
        maf_file, bam_type, alt_count_col, total_count_col = select_maf_file_and_allele_cols(sample, genotyped_mafs)
        if not maf_file or not os.path.exists(maf_file):
            print(f"[WARN] MAF file not found for sample {sample_name}, skipping {maf_file}")
            continue

        maf_df = read_maf(maf_file)
        maf_counts = get_variant_counts(maf_df, alt_count_col, total_count_col)

        # Ensure consistent types
        maf_counts["Chromosome"] = maf_counts["Chromosome"].astype(str)
        maf_counts["Start_Position"] = maf_counts["Start_Position"].astype(int)
        maf_counts["End_Position"] = maf_counts["End_Position"].astype(int)
        maf_counts["bam_type"] = bam_type
        maf_counts["tumor_normal"] = sample["tumor_normal"]
        maf_counts["assay"] = sample["assay_type"].split("_")[-1].upper()  # e.g., "access" or "impact"
        maf_counts["source"] = sample["assay_type"].split("_")[0]  # e.g., "clinical" or "research"
        maf_counts["sample_id"] = sample_name

        maf_cols=["sample_id", "tumor_normal", "assay", "source", "bam_type",
                "Chromosome","Start_Position","End_Position","Reference_Allele",
                "Tumor_Seq_Allele2","alt_count","total_count","vaf"]

        # Merge union variants with sample counts
        merged = pd.merge(
            union_df,
            maf_counts[maf_cols],
            on=["Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2"],
            how="left"
        )
        combined_list.append(merged)

    if combined_list:
        combined_df = pd.concat(combined_list, ignore_index=True)
    else:
        combined_df = union_df.copy()
        combined_df["sample_id"] = np.nan
        combined_df["patient_id"] = np.nan
        combined_df["alt_count"] = np.nan
        combined_df["total_count"] = np.nan
        combined_df["vaf"] = np.nan

    combined_df["patient_id"] = patient_id

    combined_df["call_status"] = combined_df.apply(
        lambda x: get_call_status(x, access_min_cov=access_min_cov, impact_min_cov=impact_min_cov), 
        axis=1
    )
    combined_df = combined_df.drop(columns=["Called_In"]).sort_values(by=["Chromosome", "Start_Position"])
    combined_df.to_csv(output_path, sep=",", index=False)
    print(f"[INFO] Combined variant table written to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate allele counts for all variants across all samples of a patient.")
    parser.add_argument("--patient_json", required=True)
    parser.add_argument("--genotyped_mafs", nargs="+", required=True, help="List of all genotyped MAF files")
    parser.add_argument("--union_calls_maf", required=True)
    parser.add_argument("--access_min_cov", type=int, default=100, help="Minimum coverage threshold for ACCESS samples")
    parser.add_argument("--impact_min_cov", type=int, default=50, help="Minimum coverage threshold for IMPACT samples")
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    aggregate_variants(
        args.patient_json,
        args.genotyped_mafs,
        args.union_calls_maf,
        args.output,
        access_min_cov=args.access_min_cov,
        impact_min_cov=args.impact_min_cov
    )
