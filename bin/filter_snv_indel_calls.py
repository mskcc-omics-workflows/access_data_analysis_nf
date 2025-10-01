import pandas as pd
import argparse

def filter_variants(calls_df, exclude_genes, exclude_classifications,
                    hotspot_cutoff=3, non_hotspot_cutoff=5,
                    vaf_ratio_threshold=2.0):
    """
    Add filter annotations:
    - excluded_gene: Hugo_Symbol is in exclude_genes
    - excluded_classification: Variant_Classification is in exclude_classifications
    - low_access_cov: variant not covered in any ACCESS samples (all low_coverage)
    - low_duplex_alt_count: max duplex_alt_count below cutoff depending on Hotspot column
    - low_tumor_to_normal_vaf_ratio: max tumor vaf / max normal vaf < threshold
    """
    df = calls_df.copy()
    df['filter'] = ''

    if df.empty:
        return df

    # --- Excluded genes ---
    if not set(['Hugo_Symbol']).issubset(df.columns):
        print(f"[ERROR] Missing required columns for filtering: Hugo_Symbol")
    elif exclude_genes:
        gene_mask = df['Hugo_Symbol'].isin([g for g in exclude_genes if g])
        df.loc[gene_mask, 'filter'] = 'excluded_gene'

    # --- Excluded classifications ---
    if not set(['Variant_Classification']).issubset(df.columns):
        print(f"[ERROR] Missing required columns for filtering: Variant_Classification")
    elif exclude_classifications:
        class_mask = df['Variant_Classification'].isin([c for c in exclude_classifications if c])
        df.loc[class_mask, 'filter'] = df.loc[class_mask].apply(
            lambda x: 'excluded_classification' if x['filter'] == ''
            else f"{x['filter']};excluded_classification", axis=1
        )

    # --- Low coverage check ---
    variant_keys = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"]
    required_cols = ["assay", "call_status"] + variant_keys
    if not set(required_cols).issubset(df.columns):
        print(f"[WARNING] Missing one or more required columns for coverage based filtering: {required_cols}. Skipping")
    else:
        access_df = df[df["assay"] == "ACCESS"]
        low_cov_keys = access_df.groupby(variant_keys)['call_status'].apply(
            lambda x: (x == 'low_coverage').all()
        )
        low_cov_keys = low_cov_keys[low_cov_keys].index
        df_indexed = df.set_index(variant_keys)
        low_coverage_mask = df_indexed.index.isin(low_cov_keys)
        df_indexed.loc[low_coverage_mask, 'filter'] = df_indexed.loc[low_coverage_mask].apply(
            lambda x: 'low_access_cov' if x['filter'] == ''
            else f"{x['filter']};low_access_cov", axis=1
        )
        df = df_indexed.reset_index()

    # --- Low duplex alt count filter ---
    required_duplex_cols = ["Clinical", "assay", "tumor_normal", "duplex_alt_count", "Hotspot"] + variant_keys
    if not set(required_duplex_cols).issubset(df.columns):
        print(f"[WARNING] Missing one or more required columns for duplex_alt_count filtering: {required_duplex_cols}. Skipping")
    else:
        subset = df[(df["Clinical"] != "Signed Out") & 
                    (df["assay"] == "ACCESS") & 
                    (df["tumor_normal"] == "tumor")]

        if not subset.empty:
            max_duplex = subset.groupby(variant_keys)['duplex_alt_count'].max().reset_index()
            max_duplex = max_duplex.rename(columns={"duplex_alt_count": "max_duplex_alt_count"})
            df = df.merge(max_duplex, on=variant_keys, how="left")

            low_duplex_mask = ((df["Hotspot"].astype(str).str.lower() == "yes") & (df["max_duplex_alt_count"] < hotspot_cutoff)) | \
                              ((df["Hotspot"].astype(str).str.lower() != "yes") & (df["max_duplex_alt_count"] < non_hotspot_cutoff))

            df.loc[low_duplex_mask, 'filter'] = df.loc[low_duplex_mask].apply(
                lambda x: 'low_duplex_alt_count' if x['filter'] == ''
                else f"{x['filter']};low_duplex_alt_count", axis=1
            )

            df.drop(columns=["max_duplex_alt_count"], inplace=True)

    # --- Low tumor-to-normal VAF ratio filter ---
    required_vaf_cols = ["Clinical", "tumor_normal", "vaf"] + variant_keys
    if not set(required_vaf_cols).issubset(df.columns):
        print(f"[WARNING] Missing one or more required columns for tumor-to-normal VAF ratio filtering: {required_vaf_cols}. Skipping")
    else:
        subset = df[df["Clinical"] != "Signed Out"]
        if not subset.empty:
            vaf_group = subset.groupby(variant_keys, group_keys=False).apply(
                lambda g: pd.Series({
                    "max_tumor_vaf": g.loc[g["tumor_normal"] == "tumor", "vaf"].max(skipna=True),
                    "max_normal_vaf": g.loc[g["tumor_normal"] == "normal", "vaf"].max(skipna=True)
                })
            ).reset_index()

            vaf_group["vaf_ratio"] = vaf_group.apply(
                lambda x: x["max_tumor_vaf"] / x["max_normal_vaf"] 
                          if pd.notna(x["max_normal_vaf"]) and x["max_normal_vaf"] > 0 
                          else float("inf"),
                axis=1
            )

            df = df.merge(vaf_group[[*variant_keys, "vaf_ratio"]], on=variant_keys, how="left")

            low_ratio_mask = (df["Clinical"] != "Signed Out") & (df["vaf_ratio"] < vaf_ratio_threshold)
            df.loc[low_ratio_mask, "filter"] = df.loc[low_ratio_mask].apply(
                lambda x: "low_tumor_to_normal_vaf_ratio" if x["filter"] == ""
                else f"{x['filter']};low_tumor_to_normal_vaf_ratio", axis=1
            )

            df.drop(columns=["vaf_ratio"], inplace=True)

    # change filter to PASS if it passes all the above conditions
    df.loc[df['filter'] == "", 'filter'] = "PASS"
    return df

def create_pivot_table(pass_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a pivot table with one row per unique variant (defined by required + optional columns)
    and one column per sample_id, containing alt_count/total_count(vaf).
    
    Parameters
    ----------
    pass_df : pd.DataFrame
        DataFrame of PASS-filtered variants. Must include at least:
        ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2",
         "sample_id", "alt_count", "total_count", "vaf"]

    Returns
    -------
    pivot_df : pd.DataFrame
        Pivoted table with one row per unique variant and one column per sample_id.
    """
    required_cols = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"]
    optional_cols = [
        "Hugo_Symbol", "Variant_Classification", "HGVSp",
        "HGVSp_Short", "HGVSc", "Clinical", "patient_id", "Hotspot", "CH"
    ]

    # Verify required columns exist
    missing_required = [c for c in required_cols if c not in pass_df.columns]
    if missing_required:
        raise ValueError(f"Missing required columns for variant uniqueness: {missing_required}")

    # Build final set of variant columns
    variant_cols = required_cols + [c for c in optional_cols if c in pass_df.columns]
    pass_df[variant_cols] = pass_df[variant_cols].fillna("")  # or ""

    # Check per-sample columns exist
    required_sample_cols = {"sample_id", "alt_count", "total_count", "vaf"}
    missing_sample_cols = required_sample_cols - set(pass_df.columns)
    if missing_sample_cols:
        raise ValueError(f"Missing required sample columns for pivot table: {missing_sample_cols}")

    # Build per-sample annotation string
    df = pass_df.copy()
    df["sample_annotation"] = df.apply(
        lambda x: f"{x['alt_count']}/{x['total_count']}({x['vaf']:.4f})", axis=1
    )
    # Pivot table: one row per variant, one column per sample_id
    pivot_df = df.pivot_table(
        index=variant_cols,
        columns="sample_id",
        values="sample_annotation",
        aggfunc=lambda x: ";".join(x)  # handle duplicates per variant+sample
    ).reset_index()
    print(pivot_df)
    # Reorder columns: variant columns first, then sample columns (sorted)
    pivot_df = pivot_df[
        variant_cols + sorted([col for col in pivot_df.columns if col not in variant_cols])
    ]

    return pivot_df

def main():
    parser = argparse.ArgumentParser(description="Filter SNV/Indel calls.")
    parser.add_argument("--variant_input", required=True, help="Path to variant input file.")
    parser.add_argument("--exclude_genes", default="", help="Comma-separated list of genes to exclude")
    parser.add_argument("--exclude_classifications", default="", help="Comma-separated list of variant classifications to exclude")
    parser.add_argument("--output", required=True, help="Path to save all variants with their filter column (unfiltered output).")
    parser.add_argument(
        "--output_final",
        required=True,
        help="Path to save PASS-filtered variants as a flat CSV (one row per variant-sample pair)."
    )
    parser.add_argument(
        "--output_final_table",
        required=False,
        help="Path to save PASS-filtered variants as a per-variant pivot table (one row per variant, one column per sample). "
             "If not provided, defaults to <output_final_basename>.table.csv"
    )
    parser.add_argument("--hotspot_cutoff", type=int, default=3, help="Cutoff for max duplex_alt_count of hotspot non-signed out variants (default: %(default)s).")
    parser.add_argument("--non_hotspot_cutoff", type=int, default=5, help="Cutoff for max duplex_alt_count of non-hotspot non-signed out variants (default: %(default)s).")
    parser.add_argument("--vaf_ratio_threshold", type=float, default=2.0, help="Minimum ratio threshold for max VAF of tumor samples to max VAF of normal samples for non-signed out variants (default: %(default)s).")
    args = parser.parse_args()

    exclude_genes = args.exclude_genes.split(",") if args.exclude_genes else []
    exclude_classifications = args.exclude_classifications.split(",") if args.exclude_classifications else []

    variant_df = pd.read_csv(args.variant_input, sep=",")
    if variant_df.empty:
        print("[WARN] No variants found in the input CSV.")
        filtered_df = variant_df.copy()
    else:
        filtered_df = filter_variants(
            variant_df, exclude_genes, exclude_classifications,
            hotspot_cutoff=args.hotspot_cutoff,
            non_hotspot_cutoff=args.non_hotspot_cutoff,
            vaf_ratio_threshold=args.vaf_ratio_threshold
        )
        filtered_df = filtered_df.sort_values(by=["Chromosome", "Start_Position"])

    filtered_df.to_csv(args.output, sep=",", index=False)
    print(f"[INFO] All variants with filter column saved to {args.output}")

    if 'filter' not in filtered_df.columns:
        print("[WARN] No filter column found, assuming all variants pass filters")
        filtered_df['filter'] = "PASS"

    filtered_df = filtered_df[filtered_df["filter"].isin(["", "PASS"])]
    filtered_df.to_csv(args.output_final, sep=",", index=False)
    print(f"[INFO] Filtered variants saved to {args.output_final}")

    # --- Create pivot table DataFrame and save ---
    pivot_df = create_pivot_table(filtered_df)
    table_path = args.output_final_table or args.output_final.replace(".csv", ".table.csv")
    pivot_df.to_csv(table_path, sep=",", index=False)
    print(f"[INFO] PASS-filtered pivot table saved to {table_path}")

if __name__ == "__main__":
    main()
