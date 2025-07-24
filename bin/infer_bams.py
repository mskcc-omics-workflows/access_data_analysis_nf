import sys
from collections import defaultdict
import json
import pandas as pd
import argparse
import csv
import re
import logging
import os

"""
Construct the BAM path by replacing placeholders in the template with values from the sample data.
"""

def get_bams(sample_data, template):

    if sample_data['assay_type'] == "research_access":
        sample_id = sample_data['sample_id']
        cmo_patient_id = "-".join(sample_id.split("-")[:2])
        anon_id, fl, sl = "", "", ""
    else:
        sample_id, cmo_patient_id = "", ""
        anon_id = sample_data['anon_id']
        fl, sl = anon_id[0], anon_id[1]

    bam_path = (
        template
        .replace("{sample_id}", sample_id)
        .replace("{cmo_patient_id}", cmo_patient_id)
        .replace("{anon_id}", anon_id)
        .replace("{anon_id_fl}", fl)
        .replace("{anon_id_sl}", sl)
    )

    if validate_bam(bam_path):
        return(str(os.path.realpath(bam_path)))
    else:
        return "MISSING_PATH"


def validate_bam(bam_path):
    """Validate existence of BAM file and index file."""
    if not os.path.isfile(os.path.realpath(bam_path)):
        print(f"[ERROR] BAM file not found: {bam_path}.")
        return False

    # Check for .bai index file
    bai_path_1 = bam_path + ".bai"
    bai_path_2 = bam_path.replace(".bam", ".bai")

    if os.path.isfile(bai_path_1) or os.path.isfile(bai_path_2):
        return True
    else:
        print(f"[WARNING] BAM index file (.bai) not found for: {bam_path}")
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BAM paths.")
    parser.add_argument("--sample_data", required=True)
    parser.add_argument("--template", required=True)
    args = parser.parse_args()

    get_bams(args.sample_data, args.template)
    