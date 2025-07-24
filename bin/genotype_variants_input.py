import os
import json
import pandas as pd
import argparse
import subprocess
from infer_bams import get_bams

""" 
Script to create input metadata table required for genotype_variants. Gets all the relevant bams for a set of samples in the input patient JSON.
The bam templates are passed into a helper function get_bams which replaces the template strings with the metadata from the JSON. 
"""

def build_input_table(patient_json, templates, all_calls_maf):
    """
    Main function to build a genotyping input table. Loads patient JSON, extracts BAM paths for the samples, 
    combines with patient and MAF metadata, and writes to a TSV output.
    """

    patient_data = load_patient_json(patient_json)
    combined_id = patient_data['combined_id']
    bam_paths = extract_bam_paths(patient_data, templates)
    
    bam_paths_df = pd.DataFrame(bam_paths)
    bam_paths_df['patient_id'] = combined_id
    bam_paths_df['maf'] = os.path.realpath(all_calls_maf)

    write_genotyping_table(bam_paths_df, combined_id)

def extract_bam_paths(patient_data, templates):
    """
    Builds a list of BAM paths for each sample in the patient. Fills out the fields depending on the assay_type and whether sample is tumor or normal.
    """
    bam_paths = []

    # go through each sample
    for sample_id, sample_data in patient_data["samples"].items():

        # create an entry for the sample that will store the relevant bam paths
        entry = {"sample_id": sample_id}
        assay = sample_data["assay_type"]
        tumor_normal = sample_data["tumor_normal"]

        # Research ACCESS: duplex bam and simplex bam from templates
        if assay == "research_access":
            entry["duplex_bam"] = get_bams(sample_data, templates["research_access_duplex_bam"])
            entry["simplex_bam"] = get_bams(sample_data, templates["research_access_simplex_bam"])
        # Clinical ACCESS duplex bam and simplex bam, or standard bam from templates
        elif assay == "clinical_access":
            if tumor_normal == "tumor":
                entry["duplex_bam"] = get_bams(sample_data, templates["clinical_access_duplex_bam"])
                entry["simplex_bam"] = get_bams(sample_data, templates["clinical_access_simplex_bam"])
            elif tumor_normal == "normal":
                entry["standard_bam"] = get_bams(sample_data, templates["clinical_access_standard_bam"])
        # Clinical IMPACT, standard bam from template
        elif assay == "clinical_impact":
            entry["standard_bam"] = get_bams(sample_data, templates["clinical_impact_standard_bam"])

        # add all the bam paths for one sample to the bam_paths list
        bam_paths.append(entry)
    
    # return the list of bam_paths, where each row is a sample and the columns are the different bam types
    return bam_paths

def load_patient_json(patient_json):
    """ Load the JSON file containing all metadata about the patient and their samples. """
    with open(patient_json) as json_file:
        return json.load(json_file)

def write_genotyping_table(df, patient_id):
    """ Save the genotyping input dataframe (with BAM paths + metadata) as a tsv. """
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