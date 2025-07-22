import os
import json
import pandas as pd
import argparse
import subprocess
from infer_bams import get_bams

"""
"""

def build_input_table(patient_json, templates, maf_results):
    with open(patient_json) as json_file:
        patient_data = json.load(json_file)

    bam_paths = []

    for sample_id, sample_data in patient_data["samples"].items():

        # if the sample is access research
        if sample_data['assay_type'] == "research_access":
            bam_paths.append({
            "sample_id": sample_id,
            "duplex_bam": get_bams(sample_data, templates["research_duplex_bam"]),
            "duplex_bai": get_bams(sample_data, templates["research_duplex_bai"]),
            "simplex_bam": get_bams(sample_data, templates["research_simplex_bam"]),
            "simplex_bai": get_bams(sample_data, templates["research_simplex_bai"])
            })

        # if the sample is access clinical
        if sample_data['assay_type'] == "clinical_access":
            
            if sample_data['tumor_normal'] == "tumor":
                bam_paths.append({
                "sample_id": sample_id,
                "duplex_bam": get_bams(sample_data, templates["clinical_access_duplex_bam"]),
                "duplex_bai": get_bams(sample_data, templates["clinical_duplex_bai"]),
                "simplex_bam": get_bams(sample_data, templates["clinical_simplex_bam"]),
                "simplex_bai": get_bams(sample_data, templates["clinical_simplex_bai"])
                })

            if sample_data['tumor_normal'] == "normal":
                bam_paths.append({
                "sample_id": sample_id,
                "standard_bam": get_bams(sample_data, templates["impact_standard_bam"]),
                "standard_bai": get_bams(sample_data, templates["impact_standard_bai"])
                })

        # if the sample is clinical impact
        if sample_data['assay_type'] == "clinical_impact":
            bam_paths.append({
            "sample_id": sample_id,
            "standard_bam": get_bams(sample_data, templates["impact_standard_bam"]),
            "standard_bai": get_bams(sample_data, templates["impact_standard_bai"])
            })
            
    bam_paths_df = pd.DataFrame(bam_paths)
    bam_paths_df['patient_id'] = patient_data['combined_id']
    bam_paths_df['maf'] = os.path.realpath(maf_results)

    genotyping_table_path = f"{patient_data['combined_id']}_genotyping_input.tsv"
    bam_paths_df.to_csv(genotyping_table_path, sep="\t", index=False)
    print(f"[INFO] Genotyping input saved to: {genotyping_table_path}")
    return genotyping_table_path

"""
def run_genotype_variants(genotyping_input, fasta_ref, threads):
    cmd = [
        "genotype_variants", "small_variants", "multiple-samples",
        "-i", genotyping_input,
        "-r", fasta_ref,
        "--filter-duplicate", "1",
        "-g", "/work/access/production/resources/tools/GetBaseCountsMultiSample/current/GetBaseCountsMultiSample",
        "-t", str(threads)
    ]
    print("[INFO] Running command:")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--patient_json", required=True)
    parser.add_argument("--fasta_ref", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--maf_results", required=True)

    # Research access BAMs
    parser.add_argument("--research_duplex_bam", required=True)
    parser.add_argument("--research_duplex_bai", required=True)
    parser.add_argument("--research_simplex_bam", required=True)
    parser.add_argument("--research_simplex_bai", required=True)

    # Clinical access BAMs
    parser.add_argument("--clinical_access_duplex_bam", required=True)
    parser.add_argument("--clinical_duplex_bai", required=True)
    parser.add_argument("--clinical_simplex_bam", required=True)
    parser.add_argument("--clinical_simplex_bai", required=True)

    # Clinical impact BAMs
    parser.add_argument("--impact_standard_bam", required=True)
    parser.add_argument("--impact_standard_bai", required=True)

    args = parser.parse_args()

    templates = {
        # Research access
        "research_duplex_bam": args.research_duplex_bam,
        "research_duplex_bai": args.research_duplex_bai,
        "research_simplex_bam": args.research_simplex_bam,
        "research_simplex_bai": args.research_simplex_bai,

        # Clinical access
        "clinical_access_duplex_bam": args.clinical_access_duplex_bam,
        "clinical_duplex_bai": args.clinical_duplex_bai,
        "clinical_simplex_bam": args.clinical_simplex_bam,
        "clinical_simplex_bai": args.clinical_simplex_bai,

        # Clinical impact
        "impact_standard_bam": args.impact_standard_bam,
        "impact_standard_bai": args.impact_standard_bai,

        # MAF
        "maf": args.maf_results
    }

    genotyping_input = build_input_table(args.patient_json, templates, args.maf_results)
    #run_genotype_variants(genotyping_input, args.fasta_ref, args.threads)