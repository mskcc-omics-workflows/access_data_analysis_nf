import pandas as pd
import argparse

def process_biometrics_results(genotype_file, sexmismatch_file):
    # --- Load inputs ---
    df = pd.read_csv(genotype_file)
    sex_df = pd.read_csv(sexmismatch_file)

    # Remove lines where sample is compared to itself
    df = df.loc[~df['ReferenceSample'].eq(df['QuerySample'])]

    # Group by ReferenceSampleGroup
    grouped = df.groupby('ReferenceSampleGroup')

    results = []

    for group_name, group_df in grouped:
        # --- Genotype metrics ---
        total_samples = group_df['ReferenceSample'].nunique()
        num_comparisons = group_df.shape[0]
        num_expected_matches = (group_df['Status'] == 'Expected Match').sum()
        num_unexpected_matches = (group_df['Status'] == 'Unexpected Match').sum()
        num_expected_mismatches = (group_df['Status'] == 'Expected Mismatch').sum()
        num_unexpected_mismatches = (group_df['Status'] == 'Unexpected Mismatch').sum()

        # Determine genotype QC status
        genotype_qc_status = "PASS" if (num_unexpected_mismatches == 0 and num_unexpected_matches == 0) else "FAIL"

        # Determine which ReferenceSamples had most unexpected results
        sample_max_unexpected = ""
        unexpected_mask = group_df['Status'].isin(['Unexpected Mismatch', 'Unexpected Match'])
        unexpected_counts = group_df.loc[unexpected_mask, 'ReferenceSample'].value_counts()
        if not unexpected_counts.empty:
            max_count = unexpected_counts.max()
            top_refs = unexpected_counts[unexpected_counts == max_count].index.tolist()
            sample_max_unexpected = ";".join(top_refs)

        # --- Sex mismatch QC ---
        # Filter sex mismatch table to rows where sample is in this group's ReferenceSample
        sex_subset = sex_df[sex_df['sample'].isin(group_df['ReferenceSample'].unique())]

        sex_mismatch_qc_status = "PASS"
        sample_sex_mismatch = ""

        if not sex_subset.empty:
            mismatched_samples = sex_subset.loc[sex_subset['sex_mismatch'] == True, 'sample'].tolist()
            if mismatched_samples:
                sex_mismatch_qc_status = "FAIL"
                sample_sex_mismatch = ";".join(mismatched_samples)

        results.append({
            "patient_id": group_name,
            "total_samples": total_samples,
            "num_comparisons": num_comparisons,
            "num_expected_matches": num_expected_matches,
            "num_unexpected_matches": num_unexpected_matches,
            "num_expected_mismatches": num_expected_mismatches,
            "num_unexpected_mismatches": num_unexpected_mismatches,
            "genotype_qc_status": genotype_qc_status,
            "sample_max_unexpected": sample_max_unexpected,
            "sex_mismatch_qc_status": sex_mismatch_qc_status,
            "sample_sex_mismatch": sample_sex_mismatch
        })

    return pd.DataFrame(results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize biometrics genotype comparison + sex mismatch results.")
    parser.add_argument("--genotype", required=True, help="Path to the biometrics genotype comparison CSV file.")
    parser.add_argument("--sexmismatch", required=True, help="Path to the sex mismatch output CSV file.")
    parser.add_argument("--output", required=True, help="Path to the output summary CSV file.")
    args = parser.parse_args()

    output_df = process_biometrics_results(args.genotype, args.sexmismatch)
    output_df.to_csv(args.output, index=False)

