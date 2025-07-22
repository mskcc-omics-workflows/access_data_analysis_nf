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
Script to infer bam paths given samples and a template.
"""

def get_bams(sample_data, template):

    if sample_data['assay_type'] == "research_access":
        cmo_sample_id = sample_data['sample_id']
        cmo_patient_id = "-".join(cmo_sample_id.split("-")[:2])
        anon_id, fl, sl = "", "", ""
    else:
        cmo_sample_id, cmo_patient_id = "", ""
        anon_id = sample_data['anon_id']
        fl, sl = anon_id[0], anon_id[1]
    
    bam_path = template \
            .replace("{cmo_sample_id}", cmo_sample_id) \
            .replace("{cmo_patient_id}", cmo_patient_id) \
            .replace("{anon_id}", anon_id) \
            .replace("{fl}", fl) \
            .replace("{sl}", sl)
    
    return(str(os.path.realpath(bam_path)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BAM paths.")
    parser.add_argument("--sample_data", required=True)
    parser.add_argument("--template", required=True)
    args = parser.parse_args()

    get_bams(args.sample_data, args.template)
    