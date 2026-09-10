#!/usr/bin/env python3
"""
Build the input table for ``biometrics extract``.

One row per sample with its standard BAM (taken straight from the Voyager
samplesheet) plus the metadata biometrics needs.
"""

import sys
import argparse

import pandas as pd

from patient_sheet import load_patient_data

REQUIRED_COLS = ["sample_name", "sample_group", "sample_type", "sample_sex", "sample_bam"]


def extract_rows(patient_data):
    rows = []
    for sample_id, sample_data in patient_data["samples"].items():
        standard_bam = sample_data.get("standard_bam", "")
        if not standard_bam:
            sys.stderr.write(
                f"[WARNING] Ignoring sample {sample_id}: no standard_bam in samplesheet.\n"
            )
            continue
        rows.append(
            {
                "sample_name": sample_id,
                "sample_group": patient_data["combined_id"],
                "sample_type": sample_data["tumor_normal"].title(),
                "sample_sex": patient_data["sex"],
                "sample_bam": standard_bam,
            }
        )
    return rows


def build_input_table(patient_sheet):
    patient_data = load_patient_data(patient_sheet)
    combined_id = patient_data["combined_id"]

    df = pd.DataFrame(extract_rows(patient_data))
    df = df.reindex(columns=REQUIRED_COLS, fill_value="NA")

    output_path = f"{combined_id}.biometrics_input.csv"
    df.to_csv(output_path, sep=",", index=False)
    print(f"[INFO] Biometrics input saved to: {output_path}")
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient_sheet", required=True, help="Per-patient samplesheet CSV")
    args = parser.parse_args()

    build_input_table(args.patient_sheet)
