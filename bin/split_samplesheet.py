#!/usr/bin/env python3
"""
Split the Voyager cohort samplesheet into one CSV per ``combined_id``.

Each output file keeps the original header and only the rows for that patient,
and is named ``<combined_id>.samplesheet.csv``. Downstream the pipeline fans out
one task per patient over these files.
"""

import argparse
import csv
from collections import OrderedDict


def split_samplesheet(samplesheet):
    with open(samplesheet, newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        if fieldnames is None or "combined_id" not in fieldnames:
            raise SystemExit("[ERROR] samplesheet must have a 'combined_id' column")

        rows_by_patient = OrderedDict()
        for row in reader:
            combined_id = (row.get("combined_id") or "").strip()
            if not combined_id:
                continue
            rows_by_patient.setdefault(combined_id, []).append(row)

    if not rows_by_patient:
        raise SystemExit("[ERROR] no rows with a 'combined_id' found in the samplesheet")

    for combined_id, rows in rows_by_patient.items():
        out_path = f"{combined_id}.samplesheet.csv"
        with open(out_path, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"[INFO] wrote {out_path} ({len(rows)} rows)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True, help="Path to the cohort samplesheet CSV")
    args = parser.parse_args()
    split_samplesheet(args.samplesheet)
