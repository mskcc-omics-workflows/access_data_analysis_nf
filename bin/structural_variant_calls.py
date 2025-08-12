import csv
import pandas as pd
import os
import re
import argparse
import json
import sys

def generate_sv_table(patient_json, research_access_sv_template, clinical_access_sv_file, clinical_impact_sv_file, access_structural_variant_gene_list):
    
    patient_data, cmo_id, dmp_id, combined_id = load_patient_data(patient_json)
    sv_columns = ["sample_id", "sv_type", "gene_1", "gene_2", "chr_1", "pos_1", "chr_2", "pos_2", "tumor_split_read_count", "tumor_paired_read_count", "tumor_read_count", "notes"]

    research_access_sv_calls = get_research_access_sv_calls(research_access_sv_template, patient_data, cmo_id, access_structural_variant_gene_list)
    clinical_access_sv_calls = get_clinical_access_sv_calls(clinical_access_sv_file, dmp_id, access_structural_variant_gene_list)
    clinical_impact_sv_calls = get_clinical_impact_sv_calls(clinical_impact_sv_file, dmp_id, access_structural_variant_gene_list)

    research_access_sv_calls_df = pd.DataFrame(research_access_sv_calls, columns=sv_columns)
    clinical_access_sv_calls_df = pd.DataFrame(clinical_access_sv_calls, columns=sv_columns)
    clinical_impact_sv_calls_df = pd.DataFrame(clinical_impact_sv_calls, columns=sv_columns)

    all_sv_calls_df = pd.concat([research_access_sv_calls_df, clinical_access_sv_calls_df, clinical_impact_sv_calls_df], ignore_index=True)

    all_sv_calls_df["patient_id"] = combined_id
    save_to_csv(all_sv_calls_df, combined_id)

def get_research_access_sv_calls(research_access_sv_template, patient_data, cmo_id, access_structural_variant_gene_list):

    research_access_sv_calls = []

    for sample_id, sample_data in patient_data["samples"].items():
        if sample_data['assay_type'] == "research_access" and sample_data['tumor_normal'] == "tumor":

            research_access_sv_path = infer_research_access_sv_path(research_access_sv_template, cmo_id, sample_id)

            if research_access_sv_path:
                with open(research_access_sv_path, 'r') as research_access_sv_file:
                    research_access_sv_data = csv.DictReader(research_access_sv_file, delimiter='\t')
                    for row in research_access_sv_data:
                        if row['Gene1'] in access_structural_variant_gene_list or row['Gene2'] in access_structural_variant_gene_list:
                            research_access_sv_calls.append({
                            "sample_id": row["TumorId"],
                            "sv_type": row["SV_Type"],
                            "gene_1": row["Gene1"], 
                            "gene_2": row["Gene2"],
                            "chr_1": row["Chr1"],
                            "pos_1": row["Pos1"],
                            "chr_2": row["Chr2"], 
                            "pos_2": row["Pos2"],
                            "tumor_split_read_count": row["SplitReadSupport"], 
                            "tumor_paired_read_count": row["PairEndReadSupport"], 
                            "tumor_read_count": row["TumorReadCount"], 
                            "notes": ""})
        
    return research_access_sv_calls

def infer_research_access_sv_path(research_access_sv_template, cmo_id, sample_id):
    research_access_sv_path = research_access_sv_template.replace("{cmo_patient_id}", cmo_id).replace("{sample_id}", sample_id)
    if is_valid_path(research_access_sv_path):
        return research_access_sv_path
    else:
        return False

def get_clinical_access_sv_calls(clinical_access_sv_file, dmp_id, access_structural_variant_gene_list):
    clinical_access_sv_calls = []
    if dmp_id == "":
        return clinical_access_sv_calls
    if is_valid_path(clinical_access_sv_file):
        with open(clinical_access_sv_file, 'r') as f:
            clinical_access_sv = csv.DictReader(f, delimiter="\t")
            clinical_access_sv_data = list(clinical_access_sv)
            
            for row in clinical_access_sv_data:
                if dmp_id in row["Sample_ID"]:
                    if row['Site1_Hugo_Symbol'] in access_structural_variant_gene_list or row['Site2_Hugo_Symbol'] in access_structural_variant_gene_list:
                        clinical_access_sv_calls.append({
                            "sample_id": row["Sample_ID"],
                            "sv_type": row["Class"],
                            "gene_1": row["Site1_Hugo_Symbol"],
                            "gene_2": row["Site2_Hugo_Symbol"],
                            "chr_1": row["Site1_Chromosome"],
                            "pos_1": row["Site1_Position"],
                            "chr_2": row["Site2_Chromosome"], 
                            "pos_2": row["Site2_Position"],
                            "tumor_split_read_count": row["Tumor_Split_Read_Count"], 
                            "tumor_paired_read_count": row["Tumor_Paired_End_Read_Count"], 
                            "tumor_read_count": row["Tumor_Read_Count"], 
                            "notes": row["Event_Info"]})

    return clinical_access_sv_calls

def get_clinical_impact_sv_calls(clinical_impact_sv_file, dmp_id, access_structural_variant_gene_list):
    clinical_impact_sv_calls = []
    if dmp_id == "":
        return clinical_impact_sv_calls
    if is_valid_path(clinical_impact_sv_file):
        with open(clinical_impact_sv_file, 'r') as f:
            clinical_impact_sv = csv.DictReader(f, delimiter="\t")
            clinical_impact_sv_data = list(clinical_impact_sv)
            
            for row in clinical_impact_sv_data:
                if dmp_id in row:
                    if row['Site1_Hugo_Symbol'] in access_structural_variant_gene_list or row['Site2_Hugo_Symbol'] in access_structural_variant_gene_list:
                        clinical_impact_sv_calls.append({
                            "sample_id": row["Sample_ID"],
                            "sv_type": row["Class"], 
                            "gene_1": row["Site1_Hugo_Symbol"], 
                            "gene_2": row["Site2_Hugo_Symbol"],
                            "chr_1": row["Site1_Chromosome"],
                            "pos_1": row["Site1_Position"],
                            "chr_2": row["Site2_Chromosome"], 
                            "pos_2": row["Site2_Position"],
                            "tumor_split_read_count": row["Tumor_Split_Read_Count"], 
                            "tumor_paired_read_count": row["Tumor_Paired_End_Read_Count"], 
                            "tumor_read_count": row["Tumor_Read_Count"], 
                            "notes": row["Event_Info"]})
    return clinical_impact_sv_calls

def is_valid_path(path):
    if not os.path.isfile(os.path.realpath(path)):
        print(f"[WARNING] file not found: {path}.")
        return False
    return True

def load_patient_data(patient_json):
    with open(patient_json) as json_file:
        patient_data = json.load(json_file)
        cmo_id = patient_data['cmo_id']
        dmp_id = patient_data['dmp_id']
        combined_id = patient_data['combined_id']
    return patient_data, cmo_id, dmp_id, combined_id

def save_to_csv(df, patient_id):
    df.to_csv(f'{patient_id}_SV_calls.csv', index=False)
    print(f'{patient_id}_SV_calls.csv has been created.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BAM paths.")
    parser.add_argument("--patient_json", required=True, help="Path to samples CSV file.")
    parser.add_argument("--research_access_sv_template", required=True)
    parser.add_argument("--clinical_access_sv_file", required=True)
    parser.add_argument("--clinical_impact_sv_file", required=True)
    parser.add_argument("--access_structural_variant_gene_list", required=True)
    args = parser.parse_args()

    access_structural_variant_gene_list = args.access_structural_variant_gene_list.split(",")

    generate_sv_table(args.patient_json, args.research_access_sv_template, args.clinical_access_sv_file, args.clinical_impact_sv_file, access_structural_variant_gene_list)
