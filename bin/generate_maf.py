import csv
import pandas as pd
import os
import re
import argparse
import json

def get_all_calls(patient_json, maf_template, dmp_calls_path):
    
    with open(patient_json) as json_file:
        patient_data = json.load(json_file)

    plasma_calls = []
    cmo_id = patient_data['cmo_id']
    dmp_id = patient_data['dmp_id']

    # get bam paths for research ACCESS samples
    for sample_id, sample_data in patient_data["samples"].items():

        if sample_data['assay_type'] == "research_access" and sample_data['tumor_normal'] == "tumor":

            maf_path = maf_template.replace("{cmo_patient_id}", cmo_id).replace("{cmo_sample_id}", sample_id)
            if not os.path.exists(maf_path):  
                print(f"File does not exist: {maf_path}")  
            else:
                with open(maf_path, 'r') as maf:
                    reader = csv.DictReader(maf, delimiter='\t')
                    for row in reader:
                        if row['Mutation_Status'] == 'GERMLINE':
                            continue
                        if not row['Status'] == "":
                            continue
                        
                        #plasma_calls.append(dict(row))
                        plasma_calls.append({'Hugo_Symbol': row['Hugo_Symbol'], 
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
                            'Variant_Classification': row['Variant_Classification'],
                            'HGVSp': row['HGVSp'],
                            'HGVSp_Short': row['HGVSp_Short']})
                            #'Status': row['Status']})
                           

    plasma_calls_df = pd.DataFrame(plasma_calls)

    dmp_calls = []

    if dmp_id == "NA":
        print(f'Patient {cmo_id} does not have a DMP ID. Not checking for clinical calls.')
    else:
        with open(dmp_calls_path, 'r') as dmp:
            dmp_data = (line for line in dmp if "sequenced_samples:" not in line)
            reader = csv.DictReader(dmp_data, delimiter='\t')
            for row in reader:
                if "sequenced_samples:" in row:
                    continue
                if row['Mutation_Status'] == 'GERMLINE':
                    continue
                if str(dmp_id) in row['Tumor_Sample_Barcode']:
                    #dmp_calls.append(dict(row))
                    dmp_calls.append({'Hugo_Symbol': row['Hugo_Symbol'], 
                                            'Chromosome': row["Chromosome"], 
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
                                            'Variant_Classification': row['Variant_Classification'],
                                            'HGVSp': row['HGVSp'],
                                            'HGVSp_Short': row['HGVSp_Short']})
        
    dmp_calls_df = pd.DataFrame(dmp_calls)

    all_small_calls = pd.concat([plasma_calls_df, dmp_calls_df], ignore_index=True)
    all_small_calls = all_small_calls.drop_duplicates(subset=['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Variant_Classification', 'Reference_Allele', 'Tumor_Seq_Allele2'], keep='first')
    
    all_small_calls_filtered = all_small_calls[
    (all_small_calls['Variant_Classification'] != "Silent") &
    (~all_small_calls['Hugo_Symbol'].str.contains("RP11-", na=False)) &
    (~all_small_calls['Variant_Classification'].str.contains("Intron", na=False))
    ]

    combined_id = patient_data['combined_id']

    all_small_calls_filtered.to_csv(f"{combined_id}_all_calls.maf", index=False, sep = "\t")
    print(f'{combined_id}_all_small_calls_filtered.csv has been created.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BAM paths.")
    parser.add_argument("--patient_json", required=True, help="Path to samples CSV file.")
    parser.add_argument("--maf_template", required=True)
    parser.add_argument("--dmp_calls_path", required=True)
    args = parser.parse_args()

    get_all_calls(args.patient_json, args.maf_template, args.dmp_calls_path)