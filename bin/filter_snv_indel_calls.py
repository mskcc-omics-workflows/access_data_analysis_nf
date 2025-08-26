import pandas as pd
import argparse


def filter_variants(calls_df, exclude_genes, exclude_classifications):
    """
    Add filter annotations:
    - are in the exclude genes list
    - are in exclude classifications list
    - are not covered in any of the ACCESS samples
    """

    df = calls_df.copy()
    # Initialize filter column
    df['filter'] = ''

    if df.empty:
        return df

    # Check for excluded genes
    if not set(['Hugo_Symbol']).issubset(df.columns):
        print(f"[ERROR] Missing required columns for filtering: Hugo_Symbol")
    elif exclude_genes:
        gene_mask = df['Hugo_Symbol'].isin([g for g in exclude_genes if g])
        df.loc[gene_mask, 'filter'] = 'excluded_gene'

    # Check for excluded classifications
    if not set(['Variant_Classification']).issubset(df.columns):
        print(f"[ERROR] Missing required columns for filtering: Variant_Classification")
    elif exclude_classifications:
        class_mask = df['Variant_Classification'].isin([c for c in exclude_classifications if c])
        # Combine with existing filters
        df.loc[class_mask, 'filter'] = df.loc[class_mask].apply(
            lambda x: 'excluded_classification' if x['filter'] == '' 
            else f"{x['filter']};excluded_classification", axis=1
        )

    # Check coverage in ACCESS samples
    variant_keys = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"]
    required_cols = ["assay", "call_status"] + variant_keys
    if not set(required_cols).issubset(df.columns):
        print(f"[WARNING] Missing one or more required columns for coverage based filtering: {required_cols}. Skipping")
    else:
        # For each variant, check if all ACCESS samples have 'low_coverage' in 'call_status'
        access_df = df[(df["assay"] == "ACCESS")]
        # Group by variant keys and check if all call_status == 'low_coverage'
        low_cov_keys = access_df.groupby(variant_keys)['call_status'].apply(
            lambda x: (x == 'low_coverage').all()
        )
        low_cov_keys = low_cov_keys[low_cov_keys].index  # Only variants where all ACCESS samples are low_coverage

        df_indexed = df.set_index(variant_keys)
        low_coverage_mask = df_indexed.index.isin(low_cov_keys)
        # Combine with existing filters
        df_indexed.loc[low_coverage_mask, 'filter'] = df_indexed.loc[low_coverage_mask].apply(
            lambda x: 'low_access_cov' if x['filter'] == '' 
            else f"{x['filter']};low_access_cov", axis=1
        )
        df = df_indexed.reset_index()

    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter SNV/Indel calls.")
    parser.add_argument("--variant_input", required=True, help="Path to variant input file.")
    parser.add_argument("--exclude_genes", default="", help="Comma-separated list of genes to exclude")
    parser.add_argument("--exclude_classifications", default="", 
                        help="Comma-separated list of variant classifications to exclude")
    parser.add_argument("--output", required=True, help="Path to save the variant CSV file with filter column added.")
    parser.add_argument("--output_final", required=True, help="Path to save the final variant CSV file which only keeps the PASS (filter='') variants.")
    args = parser.parse_args()
    
    # Handle empty arguments properly
    exclude_genes = args.exclude_genes.split(",") if args.exclude_genes else []
    exclude_classifications = args.exclude_classifications.split(",") if args.exclude_classifications else []
    output = args.output
    output_final = args.output_final

    # Read the input variant file
    variant_df = pd.read_csv(args.variant_input, sep=",")
    if variant_df.empty:
        print("[WARN] No variants found in the input CSV.")
        filtered_df = variant_df.copy()
    else:
        # Filter the variants based on the provided criteria
        filtered_df = filter_variants(variant_df, exclude_genes, exclude_classifications)
        filtered_df = filtered_df.sort_values(by=["Chromosome", "Start_Position"])
    # Save the filtered DataFrame to a new CSV file
    filtered_df.to_csv(output, sep=",", index=False)

    if 'filter' not in filtered_df.columns:
        print("[WARN] No filter column found, assuming all variants pass filters")
        filtered_df['filter'] = ""
    filtered_df = filtered_df[filtered_df["filter"].isin(["", "PASS"])]
    filtered_df.to_csv(output_final, sep=",", index=False)
    print(f"[INFO] All variants with filter column saved to {output}")
    print(f"[INFO] Filtered variants saved to {output_final}")
