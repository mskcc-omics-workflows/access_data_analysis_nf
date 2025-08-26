#!/usr/bin/env python3

import pandas as pd
import argparse

def read_variant_file(file_path, sep=",", column_map=None):
    """
    Read variant file and ensure required columns exist.
    Args:
        file_path (str): Path to variant file
        sep (str): File separator (default: comma)
        column_map (dict): Optional dictionary mapping file columns to required column names 
    Returns:
        pd.DataFrame: DataFrame with standardized column names
    """
    df = pd.read_csv(file_path, sep=sep)
    # Rename columns if mapping provided
    if column_map:
        df = df.rename(columns=column_map)  
    required_cols = [
        "Chromosome", "Start_Position", 
        "Reference_Allele", "Tumor_Seq_Allele2"
    ]
    # Verify required columns exist
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Available columns: {df.columns.tolist()}")
        raise ValueError(f"Missing required columns: {missing_cols}")    
    return df

def annotate_variants(variants_df, annotation_df, column_name):
    """
    Generic function to annotate variants based on overlap with a reference list.
    
    Args:
        variants_df: DataFrame containing variants to annotate
        annotation_df: DataFrame containing reference variants
        column_name: Name of the new annotation column to add
    
    Returns:
        DataFrame with new annotation column added
    """
    # Ensure consistent data types for joining
    merge_cols = [
        "Chromosome", "Start_Position", 
        "Reference_Allele", "Tumor_Seq_Allele2"
    ]
    
    for df in [variants_df, annotation_df]:
        df["Chromosome"] = df["Chromosome"].astype(str)
        df["Start_Position"] = df["Start_Position"].astype(int)
        df["Reference_Allele"] = df["Reference_Allele"].astype(str)
        df["Tumor_Seq_Allele2"] = df["Tumor_Seq_Allele2"].astype(str)

    # Check for overlap
    variants_df[column_name] = variants_df.merge(
        annotation_df[merge_cols], 
        on=merge_cols, 
        how='left',
        indicator=True
    )["_merge"].map({"left_only": "", "both": "yes"})
    
    return variants_df

def main():
    parser = argparse.ArgumentParser(description="Annotate variants with Hotspot and CH status")
    parser.add_argument("--variant_input", required=True, help="Input variant CSV file")
    parser.add_argument("--hotspot_list", required=True, help="Hotspot variants TSV file")
    parser.add_argument("--ch_list", required=True, help="CH variants TSV file")
    parser.add_argument("--output", required=True, help="Output CSV file path")
    args = parser.parse_args()

    # Read input files
    print("Reading input files...")
    variants_df = read_variant_file(args.variant_input)
    hotspot_df = read_variant_file(args.hotspot_list, sep="\t")
    ch_column_map = {
        "Chrom": "Chromosome",
        "Start": "Start_Position",
        "Ref": "Reference_Allele",
        "Alt": "Tumor_Seq_Allele2"
    }
    ch_df = read_variant_file(args.ch_list, sep="\t", column_map=ch_column_map)
    # Annotate variants with both Hotspot and CH status
    print("Annotating variants...")
    variants_df = annotate_variants(variants_df, hotspot_df, "Hotspot")
    variants_df = annotate_variants(variants_df, ch_df, "CH")

    # Write output
    print(f"Writing output to {args.output}...")
    variants_df.to_csv(args.output, index=False)
    print("Done!")

if __name__ == "__main__":
    main()
