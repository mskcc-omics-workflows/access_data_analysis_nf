import csv
import pandas as pd
import os
import re
import argparse
import json
import numpy as np

def getCNACalls(patient_json, research_access_cna_template, clinical_cna_file):

    access_cna_genes = ('AR','BRCA1','BRCA2','CDK4','EGFR','ERBB2','MET','MDM2','MLH1','MSH2','MSH6','MYC')

    with open(patient_json) as json_file:
        patient_data = json.load(json_file)

    dmp_CNA_calls = []
    cmo_id = patient_data['cmo_id']
    dmp_id = patient_data['dmp_id']

    if dmp_id == "NA":
        print(f'Patient {cmo_id} does not have a DMP ID. Not checking for clinical CNA calls.')
    else:
        with open(clinical_cna_file, 'r') as dmp_CNA:
            reader = csv.reader(dmp_CNA, delimiter="\t")
            dmp_data = list(reader)
            sample_cols = dmp_data[0]
            for index, sample_col in enumerate(sample_cols[1:], start=1):
                if str(dmp_id) in sample_col:
                    for row in dmp_data[1:]:
                        fc = row[index]
                        if fc and fc!= '0':
                            hugo_symbol = row[0]
                            if hugo_symbol in access_cna_genes:
                                dmp_CNA_calls.append({"cmo_patient_id": cmo_id, "dmp_patient_id": dmp_id, "sample_id": sample_col, "hugo_symbol": hugo_symbol, "site": "",  "fold_change": fc, "p_val": ""})
 
    if not dmp_CNA_calls:
        dmp_CNA_calls_df = pd.DataFrame(columns=["cmo_patient_id", "dmp_patient_id", "sample_id", "hugo_symbol", "site", "fold_change", "p_val", "CNA_tumor"])
    else:
        dmp_CNA_calls_df = pd.DataFrame(dmp_CNA_calls)
        dmp_CNA_calls_df["fold_change"] = pd.to_numeric(dmp_CNA_calls_df["fold_change"])
        dmp_CNA_calls_df["CNA_tumor"] = np.where(dmp_CNA_calls_df["fold_change"] > 0, "AMP", "HOMDEL")

    plasma_CNA_calls = []

    for sample_id, sample_data in patient_data["samples"].items():

        if sample_data['assay_type'] == "research_access" and sample_data['tumor_normal'] == "tumor":

            cna_path = research_access_cna_template.replace("{cmo_patient_id}", cmo_id).replace("{sample_id}", sample_id)
            if not os.path.exists(cna_path):  
                print(f"File does not exist: {cna_path}")  
            else:
                with open(cna_path, 'r') as cna:
                    reader = csv.DictReader(cna, delimiter='\t')
                    for row in reader:
                        if float(row["p.adj"]) < float(0.05):
                            sample_id = row["sample"]
                            hugo = row["region"]
                            fc = row['fc']
                            p_val = row["p.adj"]
                            site = row["Cyt"]

                            plasma_CNA_calls.append({"cmo_patient_id": cmo_id, "dmp_patient_id": dmp_id, "sample_id": sample_id, "hugo_symbol": hugo, "site": site, "fold_change": fc, "p_val": p_val})

    if not plasma_CNA_calls:
        plasma_CNA_calls_df = pd.DataFrame(columns=["cmo_patient_id", "dmp_patient_id", "sample_id", "hugo_symbol", "site", "fold_change", "p_val", "CNA_tumor"])
    else:
        plasma_CNA_calls_df = pd.DataFrame(plasma_CNA_calls)
        plasma_CNA_calls_df["CNA_tumor"] = np.where(plasma_CNA_calls_df["fold_change"] > "0", "AMP", "HOMDEL")

    all_CNA_calls = pd.concat([plasma_CNA_calls_df, dmp_CNA_calls_df], ignore_index=True)
    all_CNA_calls.to_csv(f"{cmo_id}_all_CNA_calls.csv", index=False)
    print(f'{cmo_id}_all_CNA_calls.csv has been created.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get CNA calls.")
    parser.add_argument("--patient_json", required=True, help="Path to patient JSON.")
    parser.add_argument("--research_access_cna_template", required=True)
    parser.add_argument("--clinical_cna_file", required=True)
    args = parser.parse_args()

    getCNACalls(args.patient_json, args.research_access_cna_template, args.clinical_cna_file)