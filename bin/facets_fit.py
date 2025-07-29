import glob
from pathlib import Path
import pandas as pd
import argparse
import json

def get_facets_files(facets_dir, patient_json, best_fit):
    patient_data = load_patient_data(patient_json)
    dmp_id = patient_data["dmp_id"]
    combined_id = patient_data["combined_id"]
    facets_files = []

    if best_fit:
        facets_files = find_best_facet_file(facets_dir, dmp_id)
    else:
        facets_files = find_all_facets_files(facets_dir, dmp_id)

    write_to_txt(facets_files, combined_id)

def find_all_facets_files(facets_dir, dmp_id):
    patient_path = Path(facets_dir) / dmp_id[:7]
    maf_pattern = str(patient_path / f"{dmp_id}*" / "default" / "*[0-9].ccf.maf")
    return sorted(glob.glob(maf_pattern))

def find_best_facet_file(facets_dir, dmp_id):

    #TODO: simplify and test best fit function
    
    manifest_files = find_manifest_files(facets_dir, dmp_id)
    manifest_data = []
    for manifest in manifest_files:
        manifest_path = Path(manifest)
        manifest_df = read_manifest(manifest_path)

        if manifest_df.empty:
            continue

        # Ensure required columns exist
        if not {"date_reviewed", "fit_name"}.issubset(manifest_df.columns):
            continue

        manifest_df["date_reviewed"] = pd.to_datetime(manifest_df["date_reviewed"], errors="coerce")
        manifest_df = manifest_df[manifest_df["date_reviewed"].notna()]

        # Add context to each row
        manifest_df["manifest_dir"] = manifest_path.parent
        manifest_data.append(manifest_df)

    if not manifest_data:
        return find_all_facets_files(facets_dir, dmp_id)

    manifest_data_df = pd.concat(manifest_data, ignore_index=True)

    # Try to find reviewed_best_fit with facets_qc == True
    reviewed_best = manifest_data_df[
        manifest_data_df["review_status"].str.contains("reviewed_best_fit", na=False)
        & (manifest_data_df["facets_qc"] == True)
        ]

    if not reviewed_best.empty:
        selected = reviewed_best.sort_values("date_reviewed", ascending=False).iloc[0]
    else:
        # Fall back to latest by date_reviewed
        selected = manifest_data_df.sort_values("date_reviewed", ascending=False).iloc[0]

    best_fit_folder = selected["manifest_dir"] / selected["fit_name"]
    maf_glob = str(best_fit_folder / "*[0-9].ccf.maf")
    maf_files = sorted(glob.glob(maf_glob))
    return [maf_files[0]] if maf_files else []
    
def find_manifest_files(facets_dir, dmp_id):
    base_path = Path(facets_dir) / dmp_id[:7]
    manifest_glob = str(base_path / f"{dmp_id}*" / "facets_review.manifest")
    manifest_files = sorted(glob.glob(manifest_glob))

    if not manifest_files:
        return []
    return manifest_files

def load_patient_data(patient_json):
    with open(patient_json) as json_file:
        patient_data = json.load(json_file)
        return patient_data

def write_to_txt(facets_files, patient_id):
    with open(f"{patient_id}_facets_fit.txt", "w") as output:
        if facets_files:
            for facet_file in facets_files:
                output.write(f"{facet_file}\n")
        else:
            output.write("MISSING\n")
        return f"{patient_id}_facets_fit.txt"

def read_manifest(manifest_path):
    with open(manifest_path, "r") as f:
        skip_rows = [i for i, line in enumerate(f) if line.startswith("#")]

    sep = "," if manifest_path.suffix == ".csv" else "\t"
    return pd.read_csv(manifest_path, sep=sep, skiprows=skip_rows, low_memory=False, keep_default_na=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find facets fit.")
    parser.add_argument("--facets_dir", required=True, help="Path to samples CSV file.")
    parser.add_argument("--patient_json", required=True)
    parser.add_argument("--best_fit", required=False)
    args = parser.parse_args()

    get_facets_files(args.facets_dir, args.patient_json, args.best_fit)

