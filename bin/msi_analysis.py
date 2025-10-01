import csv
import pandas as pd
import argparse
import json
import os

def parse_assay_info(sample_meta):
    """Extract source and assay from sample metadata."""
    assay_type = sample_meta.get("assay_type", "")
    parts = assay_type.split("_")
    source = parts[0] if len(parts) > 0 else ""
    assay = parts[1].upper() if len(parts) > 1 else ""
    return source, assay

def load_patient_data(patient_json):
    with open(patient_json) as json_file:
        patient_data = json.load(json_file)
        cmo_id = patient_data['cmo_id']
        dmp_id = patient_data['dmp_id']
        combined_id = patient_data['combined_id']
    return patient_data, cmo_id, dmp_id, combined_id

def infer_research_access_msi_path(research_access_msi_template, cmo_id, sample_id):
    research_access_msi_path = research_access_msi_template.replace("{cmo_patient_id}", cmo_id).replace("{sample_id}", sample_id)
    if not os.path.isfile(os.path.realpath(research_access_msi_path)):
        print(f"[WARNING] MSI file not found: {research_access_msi_path}.")
        return None
    return research_access_msi_path

def get_research_access_msi_data(patient_data, cmo_id, research_access_msi_template):
    research_access_msi_scores = []
    for sample_id, sample_data in patient_data["samples"].items():
        if sample_data.get('assay_type') == "research_access" and sample_data.get('tumor_normal') == "tumor":
            research_access_msi_path = infer_research_access_msi_path(research_access_msi_template, cmo_id, sample_id)
            if not research_access_msi_path:
                continue
            with open(research_access_msi_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    research_access_msi_scores.append({
                        'sample_id': row['Tumor_Sample_ID'],
                        'MSI_score': row['Distance_from_boundary'],
                        'coverage': row['MSI_Coverage_Tumor'],
                        'MSI_status': row['MSI_Status'],
                        'source': "research",
                        'assay': "ACCESS"
                    })
    return research_access_msi_scores

def get_clinical_access_msi_data(patient_data, dmp_id, clinical_access_msi_file):
    clinical_access_msi_scores = []
    if not dmp_id:
        return clinical_access_msi_scores
    samples = patient_data["samples"]
    # Only add samples in msi file if in the json sample list and assay_type is clinical_access and assay parsed is ACCESS
    with open(clinical_access_msi_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row.get('DMP_ASSAY_ID')
            sample_meta = samples.get(sample_id)
            if (
                sample_meta is not None and
                sample_meta.get('assay_type') == "clinical_access"
            ):
                source, assay = parse_assay_info(sample_meta)
                if assay == "ACCESS":
                    clinical_access_msi_scores.append({
                        'sample_id': sample_id,
                        'MSI_score': row.get('distance_from_boundary'),
                        'coverage': row.get('Coverage', row.get('Coverag', '')),
                        'MSI_status': row.get('msi_status'),
                        'MSI_class_type': row.get('msi_class_type', ''),
                        'note': row.get('interpretation', ''),
                        'source': "clinical",
                        'assay': "ACCESS"
                    })
    return clinical_access_msi_scores

def get_clinical_impact_msi_data(patient_data, dmp_id, clinical_impact_msi_file):
    clinical_impact_msi_scores = []
    if not dmp_id:
        return clinical_impact_msi_scores
    samples = patient_data["samples"]
    # Only add samples in msi file if in the json sample list and assay_type is clinical_impact
    with open(clinical_impact_msi_file, 'r') as f:
        filtered = (line for line in f if not line.lstrip().startswith('#'))
        reader = csv.DictReader(filtered, delimiter='\t')
        for row in reader:
            sample_id = row.get('SAMPLE_ID')
            sample_meta = samples.get(sample_id)
            if (
                sample_meta is not None and
                sample_meta.get('assay_type') == "clinical_impact"
            ):
                clinical_impact_msi_scores.append({
                    'sample_id': sample_id,
                    'MSI_score': row.get('MSI_SCORE'),
                    'coverage': row.get('SAMPLE_COVERAGE'),
                    'MSI_status': row.get('MSI_TYPE'),
                    'note': row.get('SO_COMMENTS', ''),
                    'tumor_purity': row.get('TUMOR_PURITY', ''),
                    'source': "clinical",
                    'assay': "IMPACT"
                })
    return clinical_impact_msi_scores

def generate_msi_table(patient_json, research_access_msi_template, clinical_access_msi_file, clinical_impact_msi_file):
    patient_data, cmo_id, dmp_id, combined_id = load_patient_data(patient_json)

    research_access_msi_scores = get_research_access_msi_data(patient_data, cmo_id, research_access_msi_template)
    clinical_access_msi_scores = get_clinical_access_msi_data(patient_data, dmp_id, clinical_access_msi_file)
    clinical_impact_msi_scores = get_clinical_impact_msi_data(patient_data, dmp_id, clinical_impact_msi_file)

    all_msi_scores_df = pd.concat(
        [pd.DataFrame(research_access_msi_scores),
         pd.DataFrame(clinical_access_msi_scores),
         pd.DataFrame(clinical_impact_msi_scores)],
        ignore_index=True
    )
    all_msi_scores_df["patient_id"] = combined_id
    all_msi_scores_df["cmo_patient_id"] = cmo_id
    all_msi_scores_df["dmp_patient_id"] = dmp_id

    final_df = all_msi_scores_df[['sample_id', 'patient_id', 'cmo_patient_id', 'dmp_patient_id', 'MSI_score', 'MSI_status', 'source', 'assay']].copy()
#    final_df.replace("NA", None, inplace=True)
    return final_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate MSI score table.")
    parser.add_argument("--patient_json", required=True, help="Path to patient JSON file.")
    parser.add_argument("--reseach_access_msi_template", required=True, help="Research access MSI template file path.")
    parser.add_argument("--clinical_access_msi_file", required=True, help="Clinical ACCESS MSI file path.")
    parser.add_argument("--clinical_impact_msi_file", required=True, help="Clinical IMPACT MSI file path.")
    parser.add_argument("--output", required=True, help="Output CSV file.")

    args = parser.parse_args()
    df = generate_msi_table(args.patient_json, args.reseach_access_msi_template, args.clinical_access_msi_file, args.clinical_impact_msi_file)
    df.to_csv(args.output, index=False)
    print(f'{args.output} has been created.')
