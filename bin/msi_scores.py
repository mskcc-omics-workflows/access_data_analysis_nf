import csv
import pandas as pd
import argparse
import json
import os

def generate_msi_table(patient_json, reseach_access_msi_template, clinical_access_msi_file, clinical_impact_msi_file):

    patient_data, cmo_id, dmp_id, combined_id = load_patient_data(patient_json)

    MSI_scores = []

    research_access_msi_scores = get_research_access_msi_data(patient_data, cmo_id, reseach_access_msi_template)
    clinical_access_msi_scores = get_clinical_access_msi_data(dmp_id, clinical_access_msi_file)
    clinical_impact_msi_scores = get_clinical_impact_msi_data(dmp_id, clinical_impact_msi_file)

    research_access_msi_scores_df = pd.DataFrame(research_access_msi_scores)
    clinical_access_msi_scores_df = pd.DataFrame(clinical_access_msi_scores)
    clinical_impact_msi_scores_df = pd.DataFrame(clinical_impact_msi_scores)

    all_msi_scores_df = pd.concat([research_access_msi_scores_df, clinical_access_msi_scores_df, clinical_impact_msi_scores_df], ignore_index=True)
    all_msi_scores_df["patient_id"] = combined_id

    save_to_csv(all_msi_scores_df, combined_id)

def get_research_access_msi_data(patient_data, cmo_id, reseach_access_msi_template):
    research_access_msi_scores = []
    for sample_id, sample_data in patient_data["samples"].items():
        if sample_data['assay_type'] == "research_access" and sample_data['tumor_normal'] == "tumor":
            research_access_msi_path = infer_research_access_msi_path(reseach_access_msi_template, cmo_id, sample_id)
            with open(research_access_msi_path, 'r') as reseach_access_msi_file:
                reader = csv.DictReader(reseach_access_msi_file, delimiter='\t')
                reseach_access_msi_data = list(reader)
                for row in reseach_access_msi_data:
                    research_access_msi_scores.append({'sample_id': row['Tumor_Sample_ID'], 'MSI_score': row['Distance_from_boundary'], 'coverage': row['MSI_Coverage_Tumor'], 'MSI_status': row['MSI_Status']})
    return research_access_msi_scores

def infer_research_access_msi_path(reseach_access_msi_template, cmo_id, sample_id):
    research_access_msi_path = reseach_access_msi_template.replace("{cmo_patient_id}", cmo_id).replace("{sample_id}", sample_id)
    print(research_access_msi_path)
    if not os.path.isfile(os.path.realpath(research_access_msi_path)):
        print(f"[WARNING] MSI file not found: {research_access_msi_path}.")
        return False
    return research_access_msi_path

def get_clinical_access_msi_data(dmp_id, clinical_access_msi_file):
    clinical_access_msi_scores = []
    if dmp_id == "":
        return clinical_access_msi_scores
    with open(clinical_access_msi_file, 'r') as clinical_access_msi:
        reader = csv.DictReader(clinical_access_msi)
        clinical_access_msi_data = list(reader)
        for row in clinical_access_msi_data:
            if str(dmp_id) in str(row['DMP_PATIENT_ID']):
                clinical_access_msi_scores.append({'sample_id': row['DMP_ASSAY_ID'], 'MSI_score': row['distance_from_boundary'], 'coverage': row['Coverag'], 'MSI_status': row['msi_status'], 'MSI_class_type': row['msi_class_type'], 'note': row['interpretation']})
    return clinical_access_msi_scores

def get_clinical_impact_msi_data(dmp_id, clinical_impact_msi_file):
    clinical_impact_msi_scores = []
    if dmp_id == "":
        return clinical_impact_msi_scores
    with open(clinical_impact_msi_file, 'r') as clinical_impact_msi:
        filtered_clinical_impact_msi = (line for line in clinical_impact_msi if not line.lstrip().startswith('#'))
        reader = csv.DictReader(filtered_clinical_impact_msi, delimiter='\t')
        clinical_impact_msi_data = list(reader)
        for row in clinical_impact_msi_data:
            if str(dmp_id) in str(row['PATIENT_ID']):
                clinical_impact_msi_scores.append({'sample_id': row['SAMPLE_ID'], 'MSI_score': row['MSI_SCORE'], 'coverage': row['SAMPLE_COVERAGE'], 'MSI_status': row['MSI_TYPE'], 'note': row['SO_COMMENTS'], 'tumor_purity': row['TUMOR_PURITY']})
    return clinical_impact_msi_scores

def load_patient_data(patient_json):
    with open(patient_json) as json_file:
        patient_data = json.load(json_file)
        cmo_id = patient_data['cmo_id']
        dmp_id = patient_data['dmp_id']
        combined_id = patient_data['combined_id']
    return patient_data, cmo_id, dmp_id, combined_id

def save_to_csv(df, patient_id):
    df.to_csv(f'{patient_id}_MSI_scores.csv', index=False)
    print(f'{patient_id}_MSI_scores.csv has been created.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BAM paths.")
    parser.add_argument("--patient_json", required=True, help="Path to samples CSV file.")
    parser.add_argument("--reseach_access_msi_template", required=True)
    parser.add_argument("--clinical_access_msi_file", required=True)
    parser.add_argument("--clinical_impact_msi_file", required=True)
    args = parser.parse_args()

    generate_msi_table(args.patient_json, args.reseach_access_msi_template, args.clinical_access_msi_file, args.clinical_impact_msi_file)


