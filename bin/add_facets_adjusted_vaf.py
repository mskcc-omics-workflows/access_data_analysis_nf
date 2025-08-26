import os
import pandas as pd
import numpy as np
import argparse


def read_facets_file_list(txt_path):
    try:
        df = pd.read_csv(txt_path, sep="\t")
        df = df[df["facets_path"].notnull()]
        df = df[~df["facets_path"].str.upper().eq("MISSING")]
        return df
    except Exception as e:
        print(f"[ERROR] Failed to read facets file list: {e}")
        return pd.DataFrame(columns=["facets_impact_sample", "facets_fit", "facets_path"])


def parse_facets_file(facets_impact_sample, facets_fit, facets_file, maf_cols):
    """Load a FACETS file, validate required columns, and return cleaned DataFrame or None"""

    required_cols = maf_cols + ['clonality', 'expected_alt_copies', 'tcn']
    
    try:
        facets_data = pd.read_csv(facets_file, sep="\t")
    except Exception as e:
        print(f"[ERROR] Failed to read {facets_file}: {e}")
        return None

    if not set(required_cols).issubset(facets_data.columns):
        print(f"[WARN] Skipping {facets_file}, missing required columns.")
        return None

    facets_data["Chromosome"] = facets_data["Chromosome"].astype(str)
    facets_data["Start_Position"] = facets_data["Start_Position"].astype(int)
    facets_data["End_Position"] = facets_data["End_Position"].astype(int)

    facets_data['facets_impact_sample'] = facets_impact_sample
    facets_data['facets_fit'] = facets_fit

    output_cols = maf_cols + ['facets_impact_sample', 'facets_fit', 'clonality', 'tcn', 'expected_alt_copies']

    facets_data_subset = facets_data[output_cols]

    return facets_data_subset

def calculate_adjusted_vaf(row): 
    """
    Calculate adjusted VAF based on clonality, total copy number (tcn),
    and expected alternate copies. Handles multiple possible VAF column names.
    """
    # Try different possible VAF column names
    vaf_columns = ['t_vaf', 'vaf', 't_alt_count/t_depth']
    vaf = None
    
    for col in vaf_columns:
        if col in row and pd.notnull(row[col]):
            try:
                vaf = float(row[col])
                break
            except (ValueError, TypeError):
                continue

    # Calculate VAF from alt count and depth if available
    if vaf is None and 't_alt_count' in row and 't_depth' in row:
        try:
            vaf = float(row['t_alt_count']) / float(row['t_depth'])
        except (ValueError, TypeError, ZeroDivisionError):
            return np.nan
            
    if vaf is None:
        return np.nan
        
    if row.get('clonality') == 'CLONAL':
        try:
            Tcn = float(row['tcn'])
            Talt = float(row['expected_alt_copies'])
            Ncn = 2  # Assuming diploid normal copy number
            adj_vaf = (vaf * Ncn) / (Talt + (Ncn - Tcn) * vaf)
            return adj_vaf
        except (ValueError, TypeError, ZeroDivisionError):
            return np.nan
    else:
        return np.nan

def merge_facets_info(variant_df, facets_df):
    """
    Add FACETs adjusted VAF information to variants that intersect with FACETS data.
    
    Args:
        variant_df (pd.DataFrame): DataFrame containing variant information
        facets_df (pd.DataFrame): DataFrame containing FACETS data
        
    Returns:
        pd.DataFrame: Merged dataframe with FACETS information added to matching variants
    """
    # Define all expected columns
    facets_cols = [
        'clonality', 
        'tcn', 
        'expected_alt_copies',
        'facets_impact_sample', 
        'facets_fit',
        'adjusted_vaf'
    ]

    if variant_df.empty or facets_df.empty:
        print("[WARN] Empty input dataframe(s)")
        # Create empty dataframe with all columns
        empty_df = pd.DataFrame(columns=list(variant_df.columns) + facets_cols)
        return empty_df

    # Columns to use for matching variants
    variant_cols = [
        "Chromosome", "Start_Position", "End_Position",
        "Reference_Allele", "Tumor_Seq_Allele2"
    ]

    # Ensure consistent data types for joining
    for df in [variant_df, facets_df]:
        df["Chromosome"] = df["Chromosome"].astype(str)
        df["Start_Position"] = df["Start_Position"].astype(int)
        df["End_Position"] = df["End_Position"].astype(int)
        df["Reference_Allele"] = df["Reference_Allele"].astype(str)
        df["Tumor_Seq_Allele2"] = df["Tumor_Seq_Allele2"].astype(str)

    # Merge variant and FACETS data
    merged_df = pd.merge(
        variant_df,
        facets_df[variant_cols + ['clonality', 'tcn', 'expected_alt_copies', 
                                 'facets_impact_sample', 'facets_fit']],
        on=variant_cols,
        how='inner'  # Only keep variants that exist in both dataframes
    )

    # Calculate adjusted VAF for matched variants
    if not merged_df.empty:
        merged_df['adjusted_vaf'] = merged_df.apply(calculate_adjusted_vaf, axis=1)
    
    return merged_df


def merge_all_facets_samples(variant_csv, facets_file_list):
    """
    Process all FACETS sample fits for the patient and merge results into a single dataframe.

    Args:
        variant_csv (str): Path to variant CSV file
        facets_file_list (str): Path to FACETS file list or DataFrame from read_facets_file_list
        
    Returns:
        pd.DataFrame: Combined dataframe with FACETS information for all samples
    """
    # Read variant data
    variant_df = pd.read_csv(variant_csv, sep=",")
    
    facets_file_list = read_facets_file_list(facets_file_list)
    
    if facets_file_list.empty:
        print("[WARN] No valid FACETS files found")
        return pd.DataFrame(columns=list(variant_df.columns) + [
            'clonality', 'tcn', 'expected_alt_copies',
            'facets_impact_sample', 'facets_fit', 'adjusted_vaf'
        ])

    # Process each FACETS sample
    all_results = []
    for _, row in facets_file_list.iterrows():
        facets_file = row['facets_path']
        if not os.path.exists(facets_file):
            print(f"[WARN] FACETS file not found: {facets_file}")
            continue
            
        # Read FACETS data
        try:
            facets_df = pd.read_csv(facets_file, sep="\t")
            facets_df['facets_impact_sample'] = row['facets_impact_sample']
            facets_df['facets_fit'] = row['facets_fit']
        except Exception as e:
            print(f"[ERROR] Failed to read FACETS file {facets_file}: {e}")
            continue
        
        # Merge with variants
        result = merge_facets_info(variant_df, facets_df)
        if not result.empty:
            all_results.append(result)
    
    # Combine all results
    if all_results:
        return pd.concat(all_results, ignore_index=True)
    else:
        # Return empty dataframe with all expected columns
        return pd.DataFrame(columns=list(variant_df.columns) + [
            'clonality', 'tcn', 'expected_alt_copies',
            'facets_impact_sample', 'facets_fit', 'adjusted_vaf'
        ])

def select_best_facets_impact(df):
    """
    Select the facets_impact_sample with the most non-NA adjusted VAFs in ACCESS tumor samples.
    If there is a tie, randomly select one.
    """
    if df.empty or 'facets_impact_sample' not in df.columns or 'adjusted_vaf' not in df.columns:
        return df
    if 'filter' not in df.columns:
        print("[WARN] No filter column found, assuming all variants pass filters")
        df['filter'] = ""
    df = df[df["filter"].isin(["", "PASS"]) | pd.isna(df["filter"])]
    access_mask = (df["assay"] == "ACCESS") & (df["tumor_normal"] == "tumor")
    # Count non-zero vaf per facets_impact_sample
    nonzero_counts = (
        df[access_mask & (df["vaf"].notnull()) & (df["vaf"] > 0)]
        .groupby("facets_impact_sample")
        .size()
    )
    if nonzero_counts.empty:
        return df
    max_count = nonzero_counts.max()
    top_samples = nonzero_counts[nonzero_counts == max_count].index.tolist()
    selected_sample = top_samples[0]  # Select the first one in case of a tie
    return df[df["facets_impact_sample"] == selected_sample]
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add FACETs adjusted VAF to SNV/Indel variant file.")
    parser.add_argument("--variant_csv", required=True, help="Path to the SNV/Indel variant CSV file.")
    parser.add_argument("--facets_file_list", required=True, help="Path to FACETs ccf.maf files list file.")
    parser.add_argument("--output_all_facets_samples", required=True, help="Path to save the output CSV with adjusted VAF calculations from all IMPACT samples with FACETS output.")
    parser.add_argument("--output_best_facets_sample", required=True, help="Path to save the output CSV file with FACETS adjusted VAF for IMPACT sample with most overlap with variants seen in ACCESS samples.")
    args = parser.parse_args()

    final_variant_table=merge_all_facets_samples(args.variant_csv, args.facets_file_list)
    final_variant_table.to_csv(args.output_all_facets_samples, sep=",", index=False)

    best_variant_table = select_best_facets_impact(final_variant_table)
    best_variant_table.to_csv(args.output_best_facets_sample, sep=",", index=False)

