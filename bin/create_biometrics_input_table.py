import os
import json
import pandas as pd
import argparse
import subprocess
from infer_bams import get_bams

""" 
Script to create input table for biometrics extract step. 
Gets all the relevant bam files and metadata for a set of samples in the input patient JSON.
The bam templates are passed into a helper function get_bams which replaces the template strings with the metadata from the JSON. 
"""

def build_input_table(patient_json, templates):
    """
    Main function to build a biometrics input table. Loads patient JSON, extracts standard BAM paths for the samples, 
    combines with patient metadata, and writes to a CSV output.
    """

    required_cols = ['sample_name', 'sample_group', 'sample_type', 'sample_sex', 'sample_bam']

    patient_data = load_patient_json(patient_json)
    combined_id = patient_data['combined_id']
    bam_paths = extract_bam_paths(patient_data, templates)
    
    bam_paths_df = pd.DataFrame(bam_paths)
    #bam_paths_df['patient_id'] = combined_id

    bam_paths_df = bam_paths_df.reindex(columns=required_cols, fill_value="NA")

    write_biometrics_table(bam_paths_df, combined_id)

def extract_bam_paths(patient_data, templates):
    """
    Builds a list of BAM paths for each sample in the patient. Fills out the fields depending on the assay_type and whether sample is tumor or normal.
    """
    bam_paths = []

    # go through each sample
    for sample_id, sample_data in patient_data["samples"].items():

        # create an entry for the sample that will store the relevant bam paths
        entry = {"sample_name": sample_id}
        entry["sample_group"] = patient_data["combined_id"]
        assay = sample_data["assay_type"]
        entry["sample_type"] = sample_data["tumor_normal"].title()
        entry["sample_sex"] = patient_data["sex"]
        # Research ACCESS: duplex bam and simplex bam from templates
        if assay == "research_access":
            entry["sample_bam"] = get_bams(sample_data, templates["research_access_standard_bam_template"])
        elif assay == "clinical_access":
            entry["sample_bam"] = get_bams(sample_data, templates["clinical_access_standard_bam_template"])
        elif assay == "clinical_impact":
            entry["sample_bam"] = get_bams(sample_data, templates["clinical_impact_standard_bam_template"])

        # add all the bam paths for one sample to the bam_paths list
        bam_paths.append(entry)
    
    # return the list of bam_paths, where each row is a sample and the columns are the different bam types
    return bam_paths

def load_patient_json(patient_json):
    """ Load the JSON file containing all metadata about the patient and their samples. """
    with open(patient_json) as json_file:
        return json.load(json_file)

def write_biometrics_table(df, patient_id):
    """ Save the biometrics input dataframe (with BAM paths + metadata) as a tsv. """
    output_path = f"{patient_id}.biometrics_input.csv"
    df.to_csv(output_path, sep=",", index=False)
    print(f"[INFO] Biometrics input saved to: {output_path}")
    return output_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--patient_json", required=True)

    bam_keys = [
        "research_access_standard_bam_template",
        "clinical_access_standard_bam_template",
        "clinical_impact_standard_bam_template"
    ]

    for key in bam_keys:
        parser.add_argument(f"--{key}", required=True)

    args = parser.parse_args()

    templates = {
        "research_access_standard_bam_template": args.research_access_standard_bam_template,
        "clinical_access_standard_bam_template": args.clinical_access_standard_bam_template,
        "clinical_impact_standard_bam_template": args.clinical_impact_standard_bam_template,
    }

    biometric_input = build_input_table(args.patient_json, templates)
