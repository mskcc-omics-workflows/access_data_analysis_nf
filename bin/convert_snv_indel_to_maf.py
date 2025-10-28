import os
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def convert_to_maf(input_csv, maf_columns):
    df = pd.read_csv(input_csv)

    if 'sample_id' in df.columns:
        df['Tumor_Sample_Barcode'] = df['sample_id']
    
    # Check for missing MAF columns
    missing_cols = [col for col in maf_columns if col not in df.columns]
    if missing_cols:
        print(f"[WARNING]: The following required MAF columns are missing from the CSV: {missing_cols}")
        for col in missing_cols:
            df[col] = ""

    # Convert MAF columns to expected types
    df['Start_Position'] = df['Start_Position'].astype(int)
    df['End_Position'] = df['End_Position'].astype(int)

    # Identify any extra columns that aren't part of the MAF spec
    extra_cols = [col for col in df.columns if col not in maf_columns]
    
    # Reorder columns: MAF columns first, extras after
    ordered_cols = maf_columns + extra_cols
    df = df[ordered_cols]

    csv_path = Path(input_csv)
    maf_path = Path(f"{csv_path.stem}.maf")

    return df, maf_path
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert SNV/Indel variant file to maf.")
    parser.add_argument("--variant_csvs", nargs="+", required=True, help="One or more input CSV files")
    args = parser.parse_args()

    MAF_COLUMNS = [
        'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
        'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification',
        'Tumor_Sample_Barcode', 'Match_Norm_Sample_Barcode', 'Variant_Type',
        'HGVSp', 'HGVSp_Short', 'HGVSc'
    ]  

    for input_csv in args.variant_csvs:
        print(f"Converting CSV file to MAF: {input_csv}")
        maf_df, maf_path = convert_to_maf(input_csv, MAF_COLUMNS)
        maf_df.to_csv(maf_path, sep="\t", index=False)
        print(f"[INFO] Converted MAF saved to {maf_path}")

