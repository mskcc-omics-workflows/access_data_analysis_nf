import os
import glob
import pandas as pd
import argparse
import json
from pathlib import Path

def generate_variant_table(patient_json, facets_file):

    patient_data = load_patient_data(patient_json)
    combined_id = patient_data['combined_id']
    cmo_id = patient_data['cmo_id']
    dmp_id = patient_data['dmp_id']

    maf_cols = [
        'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
        'Variant_Classification', 'Reference_Allele', 'Tumor_Seq_Allele2'
    ]

    all_variants = []

    patient_mafs = get_patient_mafs(cmo_id, dmp_id)
    facets_list = read_facets_file_list(facets_file)

    for maf in patient_mafs:

        maf_data = get_reads_from_maf(maf, maf_cols)
        if maf_data is None or maf_data.empty:
            continue
        
        if facets_list:
            for facets_file in facets_list:
                facets_data = parse_facets_file(facets_file, maf_cols)

                if facets_data is None or facets_data.empty:
                    for col in ['facets_impact_sample', 'facets_fit', 'clonality', 'tcn', 'expected_alt_copies']:
                        maf_data[col] = "NA"
                    all_variants.append(maf_data)
                    continue

                merged_data = maf_data.merge(facets_data, on=maf_cols, how='left', suffixes=('', '_facets'))
                all_variants.append(merged_data)
        else:
            for col in ['facets_impact_sample', 'facets_fit', 'clonality', 'tcn', 'expected_alt_copies']:
                maf_data[col] = "NA"
            all_variants.append(maf_data)
        

    all_variants_df = pd.concat(all_variants, ignore_index=True)

    #calculate adjusted VAF
    all_variants_df["adjusted_VAF"] = all_variants_df.apply(calculate_adjusted_vaf, axis=1)

    all_variants_df.to_csv(f"{combined_id}_SNV_table.csv", index=False)

def parse_facets_file(facets_file, maf_cols):
    """Load a FACETS file, validate required columns, and return cleaned DataFrame or None"""

    required_cols = maf_cols + ['clonality', 'expected_alt_copies', 'tcn']
    
    try:
        facets_data = pd.read_csv(facets_file, sep="\t")
    except Exception as e:
        print(f"[ERROR] Failed to read {facets_file}: {e}")
        return None, None

    if not set(required_cols).issubset(facets_data.columns):
        print(f"[WARN] Skipping {facets_file}, missing required columns.")
        return None, None

    facets_data["Chromosome"] = facets_data["Chromosome"].astype(str)
    facets_data["Start_Position"] = facets_data["Start_Position"].astype(int)
    facets_data["End_Position"] = facets_data["End_Position"].astype(int)

    facets_data['facets_impact_sample'] = Path(facets_file).parent.parent.name
    facets_data['facets_fit'] = Path(facets_file).parent.name

    output_cols = maf_cols + ['facets_impact_sample', 'facets_fit', 'clonality', 'tcn', 'expected_alt_copies']

    facets_data_subset = facets_data[output_cols]

    return facets_data_subset

def get_reads_from_maf(maf, maf_cols):
    maf_data = pd.read_csv(maf, sep='\t')

    if "DUPLEX" in maf and "SIMPLEX" in maf and "ORG" not in maf:
        sample_name = os.path.basename(maf).replace("-SIMPLEX-DUPLEX_genotyped.maf", "")
        # Compute VAF
        maf_data['t_alt_count'] = maf_data["t_alt_count_fragment_simplex_duplex"]
        maf_data['t_total_count'] = maf_data["t_total_count_fragment_simplex_duplex"]
        maf_data['VAF'] = maf_data['t_alt_count'] / maf_data['t_total_count']
        maf_data['read_type'] = "SIMPLEX-DUPLEX"

    elif "ORG-STD" in maf:
        sample_name = os.path.basename(maf).replace("_genotyped.maf", "")
        # Compute VAF
        maf_data['t_alt_count'] = maf_data["t_alt_count_standard"]
        maf_data['t_total_count'] = maf_data["t_total_count_standard"]
        maf_data['VAF'] = maf_data['t_alt_count'] / maf_data['t_total_count']
        maf_data['read_type'] = "ORG-STD"

    else:
        return None
    
    maf_data["Chromosome"] = maf_data["Chromosome"].astype(str)
    maf_data["Start_Position"] = maf_data["Start_Position"].astype(int)
    maf_data["End_Position"] = maf_data["End_Position"].astype(int)

    maf_data['sample_name'] = sample_name
    output_cols = ['sample_name'] + maf_cols + ['t_alt_count', 't_total_count', 'VAF', 'read_type']
    maf_data_subset = maf_data[output_cols]

    return maf_data_subset

def load_patient_data(patient_json):
    with open(patient_json) as json_file:
        patient_data = json.load(json_file)
    return patient_data

def write_SNV_table(SNV_table, patient_id):
    SNV_table.to_csv(f"{patient_id}_SNV_table.csv", index=False)
    print(f"{patient_id}_SNV_table.csv created")

def get_patient_mafs(cmo_id, dmp_id):
    patient_mafs = []
    for file in os.listdir("../../../results_unfilter/intermediary/genotyped_mafs"):
        if (cmo_id in file or dmp_id in file) and file.endswith(".maf"):
            patient_mafs.append((os.path.join("../../../results_unfilter/intermediary/genotyped_mafs", file)))
    return patient_mafs

def read_facets_file_list(txt_path):
    try:
        df = pd.read_csv(txt_path, sep="\t", usecols=["facets_path"])
        df = df[df["facets_path"].notnull()]
        df = df[~df["facets_path"].str.upper().eq("MISSING")]
        return df["facets_path"].tolist()
    except Exception as e:
        print(f"[ERROR] Failed to read facets file list: {e}")
        return []

def calculate_adjusted_vaf(row): 
    """
    Calculate adjusted VAF based on clonality, total copy number (tcn),
    and expected alternate copies.
    """

    try:
        vaf = float(row['VAF'])
    except (ValueError, TypeError):
        return None

    if row.get('clonality') == 'CLONAL':
        try:
            t = float(row['tcn'])
            exp = float(row['expected_alt_copies'])
            ncn = 2
            # adj_vaf = (n* vaf_value) / (exp + (n - t))
            adj_vaf = ( vaf * ncn ) / (exp + (ncn - t) * vaf )
            return adj_vaf
        except (ValueError, TypeError, ZeroDivisionError):
            return vaf

    else:
        return vaf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter calls.")
    parser.add_argument("--patient_json", required=True, help="Path to samples CSV file.")
    parser.add_argument("--facets_file", required=True, help="Path to samples CSV file.")

    args = parser.parse_args()

    generate_variant_table(args.patient_json, args.facets_file)