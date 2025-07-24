import os
import json
import pandas as pd
import argparse
import subprocess
from infer_bams import get_bams

def build_input_table(patient_json, templates, all_calls_maf):
    patient_data = load_patient_json(patient_json)
    combined_id = patient_data['combined_id']
    bam_paths = extract_bam_paths(patient_data, templates)
            
    bam_paths_df = pd.DataFrame(bam_paths)
    bam_paths_df['patient_id'] = combined_id
    bam_paths_df['maf'] = os.path.realpath(all_calls_maf)

    write_genotyping_table(bam_paths_df, combined_id)

def extract_bam_paths(patient_data, templates):
    bam_paths = []

    for sample_id, sample_data in patient_data["samples"].items():

        entry = {"sample_id": sample_id}
        assay = sample_data["assay_type"]
        tumor_normal = sample_data["tumor_normal"]

        if assay == "research_access":
            entry["duplex_bam"] = get_bams(sample_data, templates["research_access_duplex_bam"])
            entry["simplex_bam"] = get_bams(sample_data, templates["research_access_simplex_bam"])
        elif assay == "clinical_access":
            if tumor_normal == "tumor":
                entry["duplex_bam"] = get_bams(sample_data, templates["clinical_access_duplex_bam"])
                entry["simplex_bam"] = get_bams(sample_data, templates["clinical_access_simplex_bam"])
            elif tumor_normal == "normal":
                entry["standard_bam"] = get_bams(sample_data, templates["clinical_access_standard_bam"])
        elif assay == "clinical_impact":
            entry["standard_bam"] = get_bams(sample_data, templates["clinical_impact_standard_bam"])

        bam_paths.append(entry)
    
    return bam_paths

def load_patient_json(patient_json):
    with open(patient_json) as json_file:
        return json.load(json_file)

def write_genotyping_table(df, patient_id):
    output_path = f"{patient_id}_genotyping_input.tsv"
    df.to_csv(output_path, sep="\t", index=False)
    print(f"[INFO] Genotyping input saved to: {output_path}")
    return output_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--patient_json", required=True)
    parser.add_argument("--all_calls_maf", required=True)

    bam_keys = [
        "research_access_duplex_bam",
        "research_access_simplex_bam",
        "clinical_access_duplex_bam",
        "clinical_access_simplex_bam",
        "clinical_access_standard_bam",
        "clinical_impact_standard_bam"
    ]

    for key in bam_keys:
        parser.add_argument(f"--{key}", required=True)

    args = parser.parse_args()

    templates = {
        "research_access_duplex_bam": args.research_access_duplex_bam,
        "research_access_simplex_bam": args.research_access_simplex_bam,
        "clinical_access_duplex_bam": args.clinical_access_duplex_bam,
        "clinical_access_simplex_bam": args.clinical_access_simplex_bam,
        "clinical_access_standard_bam": args.clinical_access_standard_bam,
        "clinical_impact_standard_bam": args.clinical_impact_standard_bam,
    }

    genotyping_input = build_input_table(args.patient_json, templates, args.all_calls_maf)