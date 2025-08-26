import glob
from pathlib import Path
import pandas as pd
import argparse
import os
from utils import load_patient_data, save_to_tsv

def get_facets_data(facets_dir, patient_json, best_fit):

    patient_data, cmo_id, dmp_id, combined_id = load_patient_data(patient_json)

    facets_cols = ["facets_impact_sample", "facets_fit", "facets_path"]
    best_facets_fits = pd.DataFrame(columns=facets_cols)

    if not dmp_id:
        formatted_output = format_output(best_facets_fits, combined_id)
        return save_to_tsv(formatted_output, combined_id, "facets_fit", "txt")

    all_facets_sample_dirs = []
    all_facets_sample_dirs = find_all_facets_sample_dirs(facets_dir, dmp_id)

    best_fits = []

    for sample_dir in all_facets_sample_dirs:
        facets_path, fit_name = find_best_facets_fit_file(sample_dir, dmp_id)
        if facets_path is not None and fit_name is not None:
            best_fits.append({
                'facets_impact_sample': Path(sample_dir).name,
                'facets_fit': fit_name,
                'facets_path': facets_path
            })

    best_facets_fits = pd.DataFrame(best_fits, columns=facets_cols)

    formatted_output = format_output(best_facets_fits, combined_id)
    return save_to_tsv(formatted_output, combined_id, "facets_fit", "txt")

def find_all_facets_sample_dirs(facets_dir, dmp_id):
    patient_path = Path(facets_dir) / dmp_id[:7]
    dir_pattern = str(patient_path / f"{dmp_id}*")
    
    # Get all matching paths
    dmp_matches = glob.glob(dir_pattern)
    
    # Return only directories
    dirs = sorted([d for d in dmp_matches if os.path.isdir(d)])
    return dirs

def find_best_facets_fit_file(facets_dir, dmp_id):

    facets_fits = []
    fit_name = "unknown"

    facets_fits = glob.glob(os.path.join(facets_dir, "**", "*[0-9].ccf.maf"), recursive=True)

    if len(facets_fits) == 0:
        return None, None 
    
    # get the manifest file
    manifest_path = os.path.join(facets_dir, "facets_review.manifest")
    if not os.path.exists(manifest_path):
        print("MANIFEST MISSING")
        return facets_fits[0], fit_name
    
    manifest_df = read_manifest(Path(manifest_path))

    if len(facets_fits) == 1:
        facets_path = facets_fits[0]
        fit_name = Path(facets_path).parent.name
        return facets_fits[0], fit_name

    # Filter rows where facets_qc == TRUE and review_status == 'reviewed_best_fit'
    filtered = manifest_df[(manifest_df['facets_qc'] == True) & (manifest_df['review_status'] == 'reviewed_best_fit')]
    
    if filtered.empty:
        # fallback: maybe no best fits? Return None or filtered by qc only
        filtered = manifest_df[manifest_df['review_status'] == "reviewed_best_fit"]
        if filtered.empty:
            filtered = manifest_df[manifest_df['facets_qc'] == True]
            if filtered.empty:
                filtered = manifest_df[manifest_df['fit_name'] == 'default']

    filtered = filtered.sort_values('date_reviewed', ascending=False)

    # Pick the top row
    best_fit = filtered.iloc[0]
    
    # Compose path to ccf.maf file based on path + fit_name
    base_path = best_fit['path']
    fit_name = best_fit['fit_name']

    facets_fits = glob.glob(os.path.join(facets_dir, f'{fit_name}', "*[0-9].ccf.maf"), recursive=True)

    return facets_fits[0], fit_name

def format_output(facets_files, patient_id):

    if facets_files.empty:
        # Write header + MISSING row
        facets_files = pd.DataFrame([["MISSING", "MISSING", "MISSING"]], columns=["facets_impact_sample", "facets_fit", "facets_path"])

    return facets_files

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

    get_facets_data(args.facets_dir, args.patient_json, args.best_fit)

