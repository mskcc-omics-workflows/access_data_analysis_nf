import csv
import pandas as pd
import os
import re
import argparse
import json
import numpy as np

def get_all_copy_number_calls(patient_json, research_access_cna_template, clinical_cna_file, access_copy_number_gene_list, research_access_copy_number_p_value_filter):

    patient_data, cmo_id, dmp_id, combined_id = load_patient_data(patient_json)
    
    clinical_copy_number_calls = get_all_clinical_copy_number_calls(dmp_id, clinical_cna_file)
    clinical_copy_number_calls['patient_id'] = combined_id

    research_copy_number_calls = get_research_access_copy_number_calls(patient_data, cmo_id, research_access_cna_template, research_access_copy_number_p_value_filter)
    research_copy_number_calls['patient_id'] = combined_id

    all_copy_number_calls = merge_and_rescue_calls(clinical_copy_number_calls, research_copy_number_calls, access_copy_number_gene_list)

    all_copy_number_calls["cmo_patient_id"] = cmo_id
    all_copy_number_calls["dmp_patient_id"] = dmp_id
    all_copy_number_calls = all_copy_number_calls[["sample_id", "patient_id", "cmo_patient_id", "dmp_patient_id", "hugo_symbol", "site", "fold_change", "p_val", "CNA_plasma", "CNA_tumor"]]

    return save_to_csv(all_copy_number_calls, combined_id, "CNA")

def merge_and_rescue_calls(clinical_copy_number_calls, research_copy_number_calls, access_copy_number_gene_list):

    if not clinical_copy_number_calls.empty:
        clinical_subset = clinical_copy_number_calls[['hugo_symbol', 'patient_id', 'CNA_tumor']]
        all_copy_number_calls = research_copy_number_calls.merge(clinical_subset, how='left', on=["hugo_symbol", "patient_id"])
    else: 
        all_copy_number_calls = research_copy_number_calls
        all_copy_number_calls['CNA_tumor'] = np.nan

    all_copy_number_calls['CNA_plasma'] = all_copy_number_calls.apply(assign_cna_category, axis=1)
    all_copy_number_calls = all_copy_number_calls[all_copy_number_calls['CNA_plasma'] != "NA"]
    all_copy_number_calls = all_copy_number_calls.drop_duplicates(subset=['cmo_patient_id', 'sample_id', 'hugo_symbol', 'site'])

    return all_copy_number_calls

def assign_cna_category(row):

    fc = float(row.fold_change)
    if fc > 1.5 and row.hugo_symbol in access_copy_number_gene_list:
        return 'AMP'
    elif fc < -1.5 and row.hugo_symbol in access_copy_number_gene_list:
        return 'HOMDEL'
    elif fc > 1.2 and pd.notna(row['CNA_tumor']):
        return 'AMP'
    elif fc < 1.2 and pd.notna(row['CNA_tumor']):
        return 'HOMDEL'
    else:
        return "NA"
    
def get_all_clinical_copy_number_calls(dmp_id, clinical_cna_file):
    
    clinical_cna_calls = []

    if not dmp_id == "":
        with open(clinical_cna_file, 'r') as cna_file:
            reader = csv.reader(cna_file, delimiter="\t")
            dmp_data = list(reader)
            sample_cols = dmp_data[0]
            for index, sample_col in enumerate(sample_cols[1:], start=1):
                if str(dmp_id) in sample_col:
                    for row in dmp_data[1:]:
                        fc = row[index]
                        if fc and fc!= '0':
                            hugo_symbol = row[0]
                            clinical_cna_calls.append({"dmp_patient_id": dmp_id, "sample_id": sample_col, "hugo_symbol": hugo_symbol, "fold_change": fc})
 
    if not clinical_cna_calls:
        clinical_cna_calls_df = pd.DataFrame(columns=["dmp_patient_id", "sample_id", "hugo_symbol", "fold_change", "CNA_tumor"])
    else:
        clinical_cna_calls_df = pd.DataFrame(clinical_cna_calls)
        clinical_cna_calls_df["fold_change"] = pd.to_numeric(clinical_cna_calls_df["fold_change"])
        clinical_cna_calls_df["CNA_tumor"] = np.where(clinical_cna_calls_df["fold_change"] > 0, "AMP", "HOMDEL")
    
    return clinical_cna_calls_df

def get_research_access_copy_number_calls(patient_data, cmo_id, research_access_cna_template, research_access_copy_number_p_value_filter):

    research_access_cna_calls = []

    for sample_id, sample_data in patient_data["samples"].items():

        if sample_data['assay_type'] == "research_access" and sample_data['tumor_normal'] == "tumor":

            research_access_cna_path = infer_research_cna_path(research_access_cna_template, cmo_id, sample_id)
            if research_access_cna_path:
                with open(research_access_cna_path, 'r') as cna_path:
                    reader = csv.DictReader(cna_path, delimiter='\t')
                    for row in reader:
                        sample_id = row["sample"]
                        hugo = row["region"]
                        fc = row['fc']
                        p_val = row["p.adj"]
                        site = row["Cyt"]

                        research_access_cna_calls.append({"cmo_patient_id": cmo_id, "sample_id": sample_id, "hugo_symbol": hugo, "site": site, "fold_change": fc, "p_val": p_val})

    if not research_access_cna_calls:
        research_access_cna_calls_df = pd.DataFrame(columns=["cmo_patient_id", "sample_id", "hugo_symbol", "site", "fold_change", "p_val", "CNA_tumor"])
    else:
        research_access_cna_calls_df = pd.DataFrame(research_access_cna_calls)
        if research_access_copy_number_p_value_filter:
            research_access_cna_calls_df = research_access_cna_calls_df[research_access_cna_calls_df["p_val"].astype(float) < float(research_access_copy_number_p_value_filter) ]
        
    return research_access_cna_calls_df

def infer_research_cna_path(research_access_cna_template, cmo_id, sample_id):
    research_cna_path = research_access_cna_template.replace("{cmo_patient_id}", cmo_id).replace("{sample_id}", sample_id)

    if not os.path.isfile(os.path.realpath(research_cna_path)):
        print(f"[WARNING] CNA file not found: {research_cna_path}.")
        return False
    
    return research_cna_path

def load_patient_data(patient_json):
    with open(patient_json) as json_file:
        patient_data = json.load(json_file)
        cmo_id = patient_data['cmo_id']
        dmp_id = patient_data['dmp_id']
        combined_id = patient_data['combined_id']
    return patient_data, cmo_id, dmp_id, combined_id

def save_to_csv(df, patient_id, var_tag):
    output_file = f'{patient_id}_{var_tag}.csv'
    df.to_csv(output_file, index=False)
    print(f'{output_file} has been created.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get CNA calls.")
    parser.add_argument("--patient_json", required=True, help="Path to patient JSON.")
    parser.add_argument("--research_access_cna_template", required=True)
    parser.add_argument("--clinical_cna_file", required=True)
    parser.add_argument("--access_copy_number_gene_list", required=True)
    parser.add_argument("--research_access_copy_number_p_value_filter", required=True)
    args = parser.parse_args()

    access_copy_number_gene_list = args.access_copy_number_gene_list.split(",")

    get_all_copy_number_calls(args.patient_json, args.research_access_cna_template, args.clinical_cna_file, access_copy_number_gene_list, args.research_access_copy_number_p_value_filter)