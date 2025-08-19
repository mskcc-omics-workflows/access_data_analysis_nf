import pandas as pd
import argparse

def process_biometrics_genotype_comparison(biometrics_genotype_comparison_file, output_file):
    # Read the input file into a pandas DataFrame
    df = pd.read_csv(biometrics_genotype_comparison_file)
    
    # Group by ReferenceSampleGroup
    grouped = df.groupby('ReferenceSampleGroup')
    
    # Initialize an empty list to hold the results
    results = []
    
    for group_name, group_df in grouped:
        # Calculate the total number of unique samples in the ReferenceSample column for each group
        total_samples = group_df['ReferenceSample'].nunique()
        
        # Calculate the number of comparisons (rows in the table) for each group
        num_comparisons = group_df.shape[0]
        
        # Calculate the number of expected matches for each group
        num_expected_matches = group_df[group_df['Status'] == 'Expected Match'].shape[0]
        
        # Calculate the number of unexpected matches for each group
        num_unexpected_matches = group_df[group_df['Status'] == 'Unexpected Match'].shape[0]
        
        # Calculate the number of expected mismatches for each group
        num_expected_mismatches = group_df[group_df['Status'] == 'Expected Mismatch'].shape[0]
        
        # Calculate the number of unexpected mismatches for each group
        num_unexpected_mismatches = group_df[group_df['Status'] == 'Unexpected Mismatch'].shape[0]
        
        # Determine the overall status for each group
        status = 'PASS' if num_unexpected_mismatches == 0 and num_unexpected_matches == 0 else 'FAIL'
        
        # Append the results to the list
        results.append({
            'patient_id': group_name,
            'total_samples': total_samples,
            'num_comparisons': num_comparisons,
            'num_expected_matches': num_expected_matches,
            'num_unexpected_matches': num_unexpected_matches,
            'num_expected_mismatches': num_expected_mismatches,
            'num_unexpected_mismatches': num_unexpected_mismatches,
            'status': status
        })
    
    # Create a DataFrame for the output
    output_df = pd.DataFrame(results)
    
    # Write the output DataFrame to a CSV file
    output_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize biometrics genotype comparison file.")
    parser.add_argument("--input", required=True, help="Path to the biometrics genotype comparison CSV file.")
    parser.add_argument("--output", required=True, help="Path to the output summary CSV file.")
    args = parser.parse_args()
    process_biometrics_genotype_comparison(args.input, args.output)