import csv
import pandas as pd
import os
import re
import argparse
import json


def get_all_calls(patient_json, research_access_mutations_maf_template, dmp_mutations_file, exclude_genes, exclude_classifications):
    """
    Load patient data, get all mutation calls (research and clinical), merge and filter them, then write results to a file.
    """
    
    # Load in the patient JSON
    patient_data = load_patient_data(patient_json)
    combined_id = patient_data['combined_id']

    # Get the research and clinical calls from corresponding MAF files
    research_calls = get_research_access_mutations(patient_data, research_access_mutations_maf_template)
    clinical_calls = get_clinical_mutations(patient_data, dmp_mutations_file)

    # Merge and filter the calls based on the exclude gene and classification lists
    all_small_calls = merge_calls(research_calls, clinical_calls)
    all_small_calls_filtered = filter_calls(all_small_calls, exclude_genes, exclude_classifications)

    # Write the final filtered output to a maf file
    write_to_maf(all_small_calls_filtered, combined_id)

def parse_mutation_file(mutations_file, assay_type, dmp_id):
    """
    Parse maf files and extract minimum MAF columns needed for genotyping.
    assay_type: "clinical" or "research" - determines filtering logic
    dmp_id: clinical patient ID used only for clinical assay filtering
    """
    
    mutations = []
    
    # return an empty list if the maf file doesn't exist
    if not os.path.exists(mutations_file):
        print(f"File does not exist: {mutations_file}")
        return mutations
    
    with open(mutations_file, 'r') as maf:

        # for clinical mafs, exclude the metadata lines
        if assay_type == "clinical":
            maf_data = (line for line in maf if "sequenced_samples:" not in line)
        else:
            maf_data = maf

        # go through each row in the maf
        reader = csv.DictReader(maf_data, delimiter='\t')
        for row in reader:

            # exclude germline mutations from both assays
            if row['Mutation_Status'] == 'GERMLINE':
                continue

            # for research, remove variants that did not pass QC 
            if assay_type == "research" and not row['Status'] == "":
                continue
            # for clinical, only look at lines for the relevant dmp_id
            if assay_type == "clinical" and str(dmp_id) not in row['Tumor_Sample_Barcode']:
                continue 

            # extract rows that match all criteria
            try:
                variant = {
                    'Hugo_Symbol': row['Hugo_Symbol'],
                    'Chromosome': row['Chromosome'],
                    'Start_Position': int(row['Start_Position']),
                    'End_Position': int(row['End_Position']),
                    'Reference_Allele': row['Reference_Allele'],
                    'Tumor_Seq_Allele1': row['Tumor_Seq_Allele1'],
                    'Tumor_Seq_Allele2': row['Tumor_Seq_Allele2'],
                    'Tumor_Sample_Barcode': row['Tumor_Sample_Barcode'],
                    'Matched_Norm_Sample_Barcode': '',
                    't_ref_count': 0,
                    't_alt_count': 0,
                    'n_ref_count': 0,
                    'n_alt_count': 0,
                    'Variant_Classification': row.get('Variant_Classification', ''),
                    'HGVSp': row.get('HGVSp', ''),
                    'HGVSp_Short': row.get('HGVSp_Short', ''),
                }
            
                mutations.append(variant)

            # throw an error if one of the rows is missing
            except Exception as e:
                print(f"Error parsing row in {maf}: {e}")

    return mutations

def filter_calls(calls_df, exclude_genes, exclude_classifications):
    """
    Filter mutation df by excluding genes and variant classifications from exclude lists.
    exclude_genes: list of gene substrings to exclude
    exclude_classifications: list of exact variant classifications to exclude
    """

    filtered_calls = calls_df.copy()

    # Exclude rows where Hugo_Symbol contains any substring from exclude_genes
    for gene in exclude_genes:
        filtered_calls = filtered_calls[~filtered_calls['Hugo_Symbol'].str.contains(gene, na=False)]

    # Exclude rows where Variant_Classification exactly matches any exclude_classifications
    filtered_calls = filtered_calls[~filtered_calls['Variant_Classification'].isin(exclude_classifications)]

    return filtered_calls
    
def get_research_access_mutations(patient_data, research_access_mutations_maf_template):
    """
    Find every research samples in the patient_data, and parse the maf file for each sample. Collects all the calls into research_mutations list.
    """

    cmo_id = patient_data['cmo_id']
    research_mutations = []

    if not cmo_id:
        return research_mutations
    else:
        # find each research sample in the patient_data
        for sample_id, sample_data in patient_data["samples"].items():
            if sample_data['assay_type'] == "research_access" and sample_data['tumor_normal'] == "tumor":
                maf_path = research_access_mutations_maf_template.replace("{cmo_patient_id}", cmo_id).replace("{sample_id}", sample_id)
                # add the variants from the maf file into the research_mutations list
                research_mutations += parse_mutation_file(maf_path, "research", "")
    
    return research_mutations


def get_clinical_mutations(patient_data, dmp_mutations_file):
    """ Get clinical mutations from the given dmp file. """
    dmp_id = patient_data['dmp_id']
    clinical_mutations = []

    if not dmp_id:
        return clinical_mutations
    else:
         clinical_mutations = parse_mutation_file(dmp_mutations_file, "clinical", dmp_id)

    return clinical_mutations


def merge_calls(research_calls, clinical_calls):
    """ Combine research and clinical calls and remove duplicates. """

    research_calls_df = pd.DataFrame(research_calls)
    clinical_calls_df = pd.DataFrame(clinical_calls)

    all_small_calls = pd.concat([research_calls_df, clinical_calls_df], ignore_index=True)
    all_small_calls = all_small_calls.drop_duplicates(subset=['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Variant_Classification', 'Reference_Allele', 'Tumor_Seq_Allele2'], keep='first')
    
    return(all_small_calls)

def write_to_maf(calls_df, patient_id):
    """ Save call information to a tab delimited file """

    calls_df.to_csv(f"{patient_id}_all_small_calls.maf", index=False, sep = "\t")
    print(f'{patient_id}_all_small_calls.maf has been created.')

def load_patient_data(patient_json):
    with open(patient_json) as json_file:
        return json.load(json_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate all called mutations MAF.")
    parser.add_argument("--patient_json", required=True)
    parser.add_argument("--research_access_mutations_maf_template", required=True)
    parser.add_argument("--dmp_mutations_file", required=True)
    parser.add_argument("--exclude_genes")
    parser.add_argument("--exclude_classifications")
    args = parser.parse_args()

    exclude_genes = args.exclude_genes.split(",")
    exclude_classifications = args.exclude_classifications.split(",")

    get_all_calls(args.patient_json, args.research_access_mutations_maf_template, args.dmp_mutations_file, exclude_genes, exclude_classifications)

