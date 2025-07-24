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
Script that constructs BAM paths by replacing placeholders in a given template template with values from the sample data.
Supports both clinical and research paths.
Expects "sample_id", "cmo_patient_id", "anon_id", "anon_id_fl" and "anon_id_sl" as the placeholders in the templates. 
"""

def get_bams(sample_data, template):
    """ Construct a BAM path from a given template and sample-specific metadata."""

    # for research access samples, pull the sample id and cmo patient id from the metadata
    if sample_data['assay_type'] == "research_access":
        sample_id = sample_data['sample_id']
        cmo_patient_id = "-".join(sample_id.split("-")[:2])
        # the clinical placeholders are set to be empty strings
        anon_id, fl, sl = "", "", ""
    else:
        # for clinical samples, pull the anon id and the first and second letter of the anon id, set research placeholders to be empty strings
        sample_id, cmo_patient_id = "", ""
        anon_id = sample_data['anon_id']
        fl, sl = anon_id[0], anon_id[1]

    # replace placeholders with the metadata
    bam_path = (
        template
        .replace("{sample_id}", sample_id)
        .replace("{cmo_patient_id}", cmo_patient_id)
        .replace("{anon_id}", anon_id)
        .replace("{anon_id_fl}", fl)
        .replace("{anon_id_sl}", sl)
    )

    # only return the bam path if it is valid
    if validate_bam(bam_path):
        return(str(os.path.realpath(bam_path)))
    else:
        return "MISSING_PATH"


def validate_bam(bam_path):
    """
    Check if a BAM file and its index file (.bai) exist at the expected location.
    Returns True if both are present, False otherwise. 
    """

    # check if the real path of the bam path exists
    if not os.path.isfile(os.path.realpath(bam_path)):
        print(f"[ERROR] BAM file not found: {bam_path}.")
        return False

    # Check for .bam.bai or .bai index file in the same directory
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
    