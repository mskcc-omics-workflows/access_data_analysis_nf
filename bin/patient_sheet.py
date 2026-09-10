#!/usr/bin/env python3
"""
Shared loader for the Voyager-prepared per-patient samplesheet.

Voyager hands the pipeline a fully-resolved samplesheet: every BAM / MAF /
variant-file path is already populated per sample, so nothing needs to be
inferred from IDs or expanded from path templates. ``SPLIT_SAMPLESHEET`` slices
the cohort sheet into one CSV per ``combined_id``; every downstream ``bin/``
script loads that slice through :func:`load_patient_data`.

The returned structure is deliberately the same shape the pipeline previously
got from the ``INFER_SAMPLES`` patient JSON (``combined_id`` / ``cmo_id`` /
``dmp_id`` / ``sex`` / ``samples``), so the per-sample loops in the analysis
scripts stay unchanged -- only the *source* of each file path moves from a
template expansion to an explicit column.
"""

import csv

# Per-sample path columns carried by the samplesheet. Missing/blank cells load
# as "" (data type absent for that sample).
SAMPLE_PATH_COLS = (
    "duplex_bam",
    "simplex_bam",
    "unfilter_bam",
    "standard_bam",
    "maf",
    "cna_file",
    "sv_file",
    "msi_file",
)

SAMPLE_META_COLS = ("sample_id", "assay_type", "tumor_normal", "anon_id", "access_version")


def load_patient_data(sheet_csv):
    """Load a per-patient samplesheet slice into the patient-data dict.

    Returns::

        {
            "combined_id": str,
            "cmo_id": str,          # "" if the patient has no CMO id
            "dmp_id": str,          # "" if the patient has no DMP id
            "sex": "M" | "F" | "",
            "samples": {
                sample_id: {
                    "sample_id", "assay_type", "tumor_normal", "anon_id",
                    "duplex_bam", "simplex_bam", "unfilter_bam", "standard_bam",
                    "maf", "cna_file", "sv_file", "msi_file",
                },
                ...
            },
        }
    """
    patient = {
        "combined_id": "",
        "cmo_id": "",
        "dmp_id": "",
        "sex": "",
        "samples": {},
    }

    with open(sheet_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row = {k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()}

            combined_id = row.get("combined_id", "")
            if not combined_id:
                continue

            if not patient["combined_id"]:
                patient["combined_id"] = combined_id
                patient["cmo_id"] = row.get("cmo_patient_id", "") or ""
                patient["dmp_id"] = row.get("dmp_patient_id", "") or ""
                patient["sex"] = row.get("sex", "") or ""

            sample_id = row.get("sample_id", "")
            if not sample_id:
                continue

            sample = {col: (row.get(col, "") or "") for col in SAMPLE_META_COLS}
            if not sample.get("anon_id"):
                sample["anon_id"] = "NA"
            for col in SAMPLE_PATH_COLS:
                sample[col] = row.get(col, "") or ""
            patient["samples"][sample_id] = sample

    return patient
