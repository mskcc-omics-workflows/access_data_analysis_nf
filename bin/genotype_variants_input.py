import os
import pandas as pd
import argparse
from infer_bams import get_bams
from utils import load_patient_data, save_to_tsv

""" 
Script to create input metadata table required for genotype_variants. Gets all the relevant bams for a set of samples in the input patient JSON.
The bam templates are passed into a helper function get_bams which replaces the template strings with the metadata from the JSON. 
"""

def build_input_table(patient_json, templates, all_calls_maf):
    """
    Main function to build a genotyping input table. Loads patient JSON, extracts BAM paths for the samples, 
    combines with patient and MAF metadata, and writes to a TSV output.
    """

    required_cols = ['patient_id', 'sample_id', 'standard_bam', 'duplex_bam', 'simplex_bam', 'maf']

    patient_data, cmo_id, dmp_id, combined_id = load_patient_data(patient_json)
    bam_paths = extract_bam_paths(patient_data, templates)
    
    bam_paths_df = pd.DataFrame(bam_paths)
    bam_paths_df['maf'] = os.path.realpath(all_calls_maf)

    bam_paths_df = bam_paths_df.reindex(columns=required_cols, fill_value="NA")

    save_to_tsv(bam_paths_df, combined_id, "genotyping_input", "tsv")

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
            if tumor_normal == "tumor":
                entry["duplex_bam"] = get_bams(sample_data, templates["research_access_duplex_bam_template"])
                entry["simplex_bam"] = get_bams(sample_data, templates["research_access_simplex_bam_template"])
            elif tumor_normal == "normal":
                entry["standard_bam"] = get_bams(sample_data, templates["research_access_unfilter_bam_template"])
        # Clinical ACCESS duplex bam and simplex bam, or standard bam from templates
        elif assay == "clinical_access":
            if tumor_normal == "tumor":
                entry["duplex_bam"] = get_bams(sample_data, templates["clinical_access_duplex_bam_template"])
                entry["simplex_bam"] = get_bams(sample_data, templates["clinical_access_simplex_bam_template"])
            elif tumor_normal == "normal":
                entry["standard_bam"] = get_bams(sample_data, templates["clinical_access_unfilter_bam_template"])
        # Clinical IMPACT, standard bam from template
        elif assay == "clinical_impact":
            entry["standard_bam"] = get_bams(sample_data, templates["clinical_impact_standard_bam_template"])

        # add all the bam paths for one sample to the bam_paths list
        bam_paths.append(entry)
    
    # return the list of bam_paths, where each row is a sample and the columns are the different bam types
    return bam_paths

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--patient_json", required=True)
    parser.add_argument("--all_calls_maf", required=True)

    bam_keys = [
        "research_access_duplex_bam_template",
        "research_access_simplex_bam_template",
        "research_access_unfilter_bam_template",
        "clinical_access_duplex_bam_template",
        "clinical_access_simplex_bam_template",
        "clinical_access_unfilter_bam_template",
        "clinical_impact_standard_bam_template"
    ]

    for key in bam_keys:
        parser.add_argument(f"--{key}", required=True)

    args = parser.parse_args()

    templates = {
        "research_access_duplex_bam_template": args.research_access_duplex_bam_template,
        "research_access_simplex_bam_template": args.research_access_simplex_bam_template,
        "research_access_unfilter_bam_template": args.research_access_unfilter_bam_template,
        "clinical_access_duplex_bam_template": args.clinical_access_duplex_bam_template,
        "clinical_access_simplex_bam_template": args.clinical_access_simplex_bam_template,
        "clinical_access_unfilter_bam_template": args.clinical_access_unfilter_bam_template,
        "clinical_impact_standard_bam_template": args.clinical_impact_standard_bam_template,
    }

    genotyping_input = build_input_table(args.patient_json, templates, args.all_calls_maf)