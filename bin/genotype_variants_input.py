#!/usr/bin/env python3
"""
Build the metadata TSV consumed by ``genotype_variants small_variants
multiple-samples``.

Every sample of the patient is genotyped against the patient's union call set
(``--all_calls_maf``). BAM paths come straight from the Voyager samplesheet:

* research / clinical ACCESS tumor  -> duplex_bam + simplex_bam
* research / clinical ACCESS normal -> standard_bam := unfilter_bam
* clinical IMPACT                   -> standard_bam
"""

import os
import sys
import argparse

import pandas as pd

from patient_sheet import load_patient_data

# genotype_variants multiple-samples reads these columns by name.
REQUIRED_COLS = ["patient_id", "sample_id", "standard_bam", "duplex_bam", "simplex_bam", "maf"]


def build_bam_entry(sample_data):
    """Return the {standard,duplex,simplex}_bam dict for one sample, or None to skip."""
    assay = sample_data["assay_type"]
    tumor_normal = sample_data["tumor_normal"]
    entry = {}

    if assay in ("research_access", "clinical_access"):
        if tumor_normal == "tumor":
            entry["duplex_bam"] = sample_data.get("duplex_bam", "")
            entry["simplex_bam"] = sample_data.get("simplex_bam", "")
        else:
            entry["standard_bam"] = sample_data.get("unfilter_bam", "")
    elif assay == "clinical_impact":
        entry["standard_bam"] = sample_data.get("standard_bam", "")
    else:
        return None

    # Drop the sample if any BAM it needs is absent from the samplesheet.
    for key, value in entry.items():
        if not value:
            sys.stderr.write(
                f"[WARNING] Ignoring sample {sample_data['sample_id']}: no {key} in samplesheet.\n"
            )
            return None
    return {k: os.path.realpath(v) for k, v in entry.items()}


def build_input_table(patient_sheet, all_calls_maf):
    patient_data = load_patient_data(patient_sheet)
    combined_id = patient_data["combined_id"]

    rows = []
    for sample_id, sample_data in patient_data["samples"].items():
        entry = build_bam_entry(sample_data)
        if entry is None:
            continue
        entry["sample_id"] = sample_id
        rows.append(entry)

    df = pd.DataFrame(rows)
    df["patient_id"] = combined_id
    df["maf"] = os.path.realpath(all_calls_maf)
    df = df.reindex(columns=REQUIRED_COLS, fill_value="")

    output_path = f"{combined_id}_genotyping_input.tsv"
    df.to_csv(output_path, sep="\t", index=False)
    print(f"[INFO] Genotyping input saved to: {output_path}")
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient_sheet", required=True, help="Per-patient samplesheet CSV")
    parser.add_argument("--all_calls_maf", required=True, help="Patient union call set MAF")
    args = parser.parse_args()

    build_input_table(args.patient_sheet, args.all_calls_maf)
