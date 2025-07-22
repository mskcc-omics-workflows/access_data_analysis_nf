import os
import glob
import pandas as pd
import argparse
import json

def filter_calls(patient_json, facet_file):

    with open(patient_json) as json_file:
        patient_data = json.load(json_file)

    facets_data = pd.read_csv(facet_file, sep="\t")
    print(facets_data.columns.tolist())


    patient_id = patient_data['cmo_patient_id']

    maf_cols = [
        'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
        'Variant_Classification', 'Reference_Allele', 'Tumor_Seq_Allele2'
    ]

    patient_mafs = []

    for file in os.listdir("../../../results/genotyping"):
        if patient_id in file and file.endswith(".maf"):
            patient_mafs.append((os.path.join("../../../results/genotyping", file)))

    
    SNV_table = pd.DataFrame(columns=maf_cols)


    for maf in patient_mafs:
        print(maf)
        maf_data = pd.read_csv(maf, sep='\t')


        if "DUPLEX" in maf and "SIMPLEX" in maf and "ORG" not in maf:

            if SNV_table.empty:
                SNV_table = pd.concat([SNV_table, maf_data[maf_cols]])
                SNV_table["Start_Position"] = SNV_table["Start_Position"].astype(int)
                SNV_table["End_Position"] = SNV_table["End_Position"].astype(int)
            # for each row in the maf
            # create a column for the sample in the SNV_table
            # calculate the VAF (t_alt_count_fragment_simplex_duplex/t_total_count_fragment_simplex_duplex)
            # use the identifying variant columns (hugo symbol, chromosome, start position, end position)
                # add the VAF to the row that corresponds to the correct variant into the new column for this sample

            sample_name = os.path.basename(maf).replace("-SIMPLEX-DUPLEX_genotyped.maf", "")

            # Compute VAF
            maf_data["reads"] = "(" + maf_data["t_alt_count_fragment_simplex_duplex"].astype(str) + "/" + maf_data["t_total_count_fragment_simplex_duplex"].astype(str) + ")"

            # Prepare mini table with variant info + VAF
            sample_vaf_df = maf_data[maf_cols + ["reads"]].copy()
            sample_vaf_df.rename(columns={"reads": sample_name}, inplace=True)
            print(sample_vaf_df)

            # Merge VAFs into SNV_table
            SNV_table = SNV_table.merge(sample_vaf_df, on=maf_cols, how='left')
            print(SNV_table)

        
        elif "-IM" in maf and "ORG-STD" in maf:
            sample_name = os.path.basename(maf).replace("_genotyped.maf", "")

            # Compute VAF
            maf_data["reads"] = "(" + maf_data["t_alt_count_standard"].astype(str) + "/" + maf_data["t_total_count_standard"].astype(str) + ")"

            # Prepare mini table with variant info + VAF
            sample_vaf_df = maf_data[maf_cols + ["reads"]].copy()
            sample_vaf_df.rename(columns={"reads": sample_name}, inplace=True)
            print(sample_vaf_df)

            # Merge VAFs into SNV_table
            SNV_table = SNV_table.merge(sample_vaf_df, on=maf_cols, how='left')
            print(SNV_table)

        else:
            continue

    SNV_table = merge_facets_info(SNV_table, facet_file)
    SNV_table['adjusted_VAF'] = SNV_table.apply(calculate_adjusted_vaf, axis=1)

    SNV_table.to_csv(f"{patient_id}_SNV_table.csv", index=False)
    print(f"{patient_id}_SNV_table.csv created")

def merge_facets_info(snv_df, facet_file):
    """Merge clonality and CN info from FACETS into SNV table"""

    print(f"[INFO] Loading facets file: {facet_file}")
    facet_df = pd.read_csv(facet_file, sep="\t")

    maf_cols = [
        'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
        'Variant_Classification', 'Reference_Allele', 'Tumor_Seq_Allele2'
    ]

    # Columns to add from facets
    extra_cols = ['clonality', 'expected_alt_copies', 'tcn']
    required_cols = maf_cols + extra_cols

    missing_cols = [col for col in required_cols if col not in facet_df.columns]
    if missing_cols:
        print(f"[WARN] Missing columns in FACETS file: {missing_cols}. Skipping merge.")
        return snv_df

    # Convert types and merge
    facet_df["Start_Position"] = facet_df["Start_Position"].astype(int)
    facet_df["End_Position"] = facet_df["End_Position"].astype(int)
    facet_subset = facet_df[required_cols].drop_duplicates()

    merged = snv_df.merge(facet_subset, on=maf_cols, how="left")
    print("[INFO] Merged SNV table with FACETS info.")
    return merged

def calculate_adjusted_vaf(row):
    """
    Calculate adjusted VAF based on clonality, total copy number (tcn),
    and expected alternate copies.
    """
    try:
        vaf_str = row['reads']  # The "reads" column has format like "(alt_count/total_count)"
        # Extract numeric VAF fraction from reads string
        if pd.isna(vaf_str) or vaf_str == '':
            return None
        parts = vaf_str.strip("()").split("/")
        alt_count = float(parts[0])
        total_count = float(parts[1])
        if total_count == 0:
            return 0.0
        vaf = alt_count / total_count
    except Exception:
        # If reads format or conversion fails, fallback
        return None

    clonality = row.get('clonality', None)
    t = row.get('tcn', None)
    exp = row.get('expected_alt_copies', None)
    ncn = 2  # normal copy number assumption

    if pd.isna(vaf) or vaf is None:
        return None

    if clonality == 'CLONAL' and pd.notna(t) and pd.notna(exp):
        try:
            t = float(t)
            exp = float(exp)
            adj_vaf = (vaf * ncn) / (exp + (ncn - t) * vaf)
            return adj_vaf
        except Exception:
            return vaf
    else:
        # For SUBCLONAL or unknown clonality, you can decide how to treat — here just return vaf
        return vaf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter calls.")
    parser.add_argument("--patient_json", required=True, help="Path to samples CSV file.")
    parser.add_argument("--facet_file", required=True, help="Path to samples CSV file.")

    args = parser.parse_args()

    filter_calls(args.patient_json, args.facet_file)