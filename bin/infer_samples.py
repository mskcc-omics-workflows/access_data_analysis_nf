import sys
import json
import argparse
import csv
import re
import os
import glob

"""
Script to read in the CMO/DMP sample IDs and retrieve all associated
research access, clinical access, and clinical impact samples.
Samples specified in the include and exclude files will be filtered as well.
Output is one JSON file per patient, containing all samples relevant to the patient under a combined cmo/dmp id.
"""

def get_all_samples(id_mapping_file, research_access_bam_dir_template, clinical_access_key_file, clinical_impact_key_file, keep_research_samples_file, exclude_samples_file, clinical_access_sample_regex_pattern, clinical_impact_sample_regex_pattern, research_access_manifest_file_template, research_access_mutations_maf_template):
    """ Main logic function to get all samples from the id mapping file, split them by patient, get the relevant samples, and save to JSON. """

    # Extract the cmo ids, dmp ids, sex, and combined ids from the input file.
    id_list = get_id_mapping(id_mapping_file)

    # Go through each patient one by one
    for patient_data in id_list:
        # Get all the ids of the patient
        combined_id = patient_data["combined_id"]
        cmo_id = patient_data["cmo_id"]
        dmp_id = patient_data["dmp_id"]
        sex = patient_data["sex"]
        sample_dict = { combined_id: patient_data }

        # 1. Find research samples if patient has a cmo id
        if cmo_id:
            # Get research samples from directory structure
            research_samples = find_research_samples(research_access_bam_dir_template, cmo_id)
            
            # 2. If a real include file is provided, filter research samples to only keep those in the include list
            # NOTE: Due to Nextflow's limitation with handling null path parameters, we use a placeholder file
            # named 'NO_INCLUDE_FILE' when no include file is actually provided. We need to detect this special
            # case and skip filtering when this placeholder is used. This is a temporary workaround that should
            # be replaced with a more elegant solution in future pipeline versions.
            if keep_research_samples_file: # and not os.path.basename(keep_research_samples_file) == "NO_INCLUDE_FILE":
                include_list = get_include_list(keep_research_samples_file, cmo_id)
                research_samples = [s for s in research_samples if s in include_list]
            
            # Add research samples to sample dictionary
            for sample_id in research_samples:
                sample_dict[combined_id]["samples"][sample_id] = {
                    "sample_id": sample_id,
                    "tumor_normal": infer_tumor_normal(sample_id),
                    "assay_type": "research_access",
                    "anon_id": "NA",
                    "access_version": infer_access_version(sample_id, research_access_manifest_file_template)
                }

                donor_id = infer_research_donor(cmo_id, sample_id, research_access_mutations_maf_template)
                if donor_id:
                    sample_dict[combined_id]["samples"][sample_id]["donor_id"] = donor_id

        # 3. Find clinical samples if patient has a dmp id
        if dmp_id:
            find_clinical_samples(clinical_impact_key_file, clinical_access_key_file, dmp_id, combined_id, 
                                sample_dict, clinical_access_sample_regex_pattern, clinical_impact_sample_regex_pattern)
        
        # 4. Exclude samples if exclude file is provided
        if exclude_samples_file:
            exclude_samples(sample_dict, exclude_samples_file, cmo_id, dmp_id, combined_id)
        
        # 5. Save the final dictionary to json
        save_to_json(sample_dict, combined_id)
        
    if not id_list:
        print("No samples found in input file.")

def find_research_samples(research_access_bam_dir_template, cmo_id):
    """ Find all valid research samples for a patient in the directory structure. """
    research_access_bam_dir_root = research_access_bam_dir_template.split("/{sample_id}")[0].replace("{cmo_patient_id}", cmo_id)
    
    research_samples = []
    try:
        for sample_name in os.listdir(research_access_bam_dir_root):
            current_path = os.path.join(research_access_bam_dir_root, sample_name, "current")
            
            # Check for valid sample (has current dir with bam files)
            if os.path.isdir(current_path) and any(f.endswith(".bam") for f in os.listdir(current_path)):
                research_samples.append(sample_name)
            else:
                print(f'{sample_name} is not a valid research sample.')
    except FileNotFoundError:
        print(f"Directory not found: {research_access_bam_dir_root}")
    
    return research_samples

def find_clinical_samples(clinical_impact_key_file, clinical_access_key_file, dmp_id, combined_id, 
                         sample_dict, clinical_access_sample_regex_pattern, clinical_impact_sample_regex_pattern):
    """ Find all clinical access and impact samples for a patient. """
    # Search the clinical key files for samples matching the dmp_id
    search_dmp_key_file(clinical_access_key_file, "clinical_access", clinical_access_sample_regex_pattern, 
                       combined_id, dmp_id, sample_dict)
    search_dmp_key_file(clinical_impact_key_file, "clinical_impact", clinical_impact_sample_regex_pattern, 
                       combined_id, dmp_id, sample_dict)

def search_dmp_key_file(dmp_key_path, assay_type, regex_pattern, combined_id, dmp_id, sample_dict):
    """ Search the key file for samples matching the dmp_id and regex pattern. """
    with open(dmp_key_path, 'r') as key_file:
        for line in key_file:
            if re.search(rf"{dmp_id}-{regex_pattern}", line):
                sample_id, anon_id, sample_type = parse_key_line(line)
                sample_dict[combined_id]["samples"][sample_id] = {
                    "sample_id": sample_id,
                    "tumor_normal": sample_type,
                    "assay_type": assay_type,
                    "anon_id": anon_id
                }

def exclude_samples(sample_dict, exclude_samples_file, cmo_id, dmp_id, combined_id):
    """ Remove samples that are in the exclude list. """
    exclude_list = []
    
    # Get all exclude samples for this patient
    if cmo_id:
        exclude_list.extend(get_exclude_list(exclude_samples_file, cmo_id))
    if dmp_id:
        exclude_list.extend(get_exclude_list(exclude_samples_file, dmp_id))
    
    # Remove samples in exclude list
    samples_to_remove = []
    for sample_id in sample_dict[combined_id]["samples"]:
        if sample_id in exclude_list:
            samples_to_remove.append(sample_id)
            print(f"Excluded sample: {sample_id}")
    
    for sample_id in samples_to_remove:
        del sample_dict[combined_id]["samples"][sample_id]

def parse_key_line(line):
    """ Parse a line from a key file to extract sample_id, anon_id, and infer sample type. """
    cols = line.split(sep=',')
    sample_id = cols[0]
    anon_id = cols[1]

    # Use the sample id to infer if sample is tumor or normal 
    sample_type = infer_tumor_normal(sample_id)

    # Remove anon id suffix from access samples
    if "-XS" in sample_id:
        anon_id = re.sub(r'-standard$', '', anon_id)

    return sample_id, anon_id, sample_type

def infer_tumor_normal(sample_id):
    """ Determine whether sample is tumor or normal, assuming sample format of C-XXXX-NXXXX or P-XXXX-NXXXX """
    parts = sample_id.split("-")
    if len(parts) >= 3 and parts[2].startswith("N"):
        return "normal"
    else:
        return "tumor"

def infer_access_version(sample_id, research_access_manifest_file_template):
    
    manifest_list = manifest_files = glob.glob(research_access_manifest_file_template)
    for manifest in manifest_list:
        with open(manifest) as f:
            for line in f:
                if sample_id in line:
                    if "MSK-ACCESS-v1" in line:
                        return "XS1"
                    if "MSK-ACCESS-v2" in line:
                        return "XS2"
                    return "not_found"

def infer_research_donor(cmo_id, sample_id, research_access_mutations_maf_template):
    
    if infer_tumor_normal(sample_id) == "normal":
        return None

    maf_pattern = research_access_mutations_maf_template.replace("{cmo_patient_id}", cmo_id).replace("{sample_id}", sample_id).replace("{donor_id}", "*")
    maf_file = glob.glob(maf_pattern)
    donor_id = maf_file[0].split('.')[1]

    return donor_id

def get_id_mapping(id_mapping_file):
    """ Read in the cmo and dmp ids, generate the combined patient id, and save each patient id in a list. """
    id_list = []

    with open(id_mapping_file, newline='') as ids:
        reader = csv.DictReader(ids)
        for row in reader:
            cmo_id = row['cmo_patient_id'].strip() if 'cmo_patient_id' in row else ""
            dmp_id = row['dmp_patient_id'].strip() if 'dmp_patient_id' in row else ""
            sex = row['sex'].strip() if 'sex' in row else ""
            combined_id = get_combined_patient_id(cmo_id, dmp_id)
            if cmo_id or dmp_id:
                id_list.append({"combined_id": combined_id, "cmo_id": cmo_id, "dmp_id": dmp_id, "sex": sex, "samples": {}})

    return id_list

def get_include_list(keep_research_samples_file, patient_id):
    """ Read in the include csv file and return samples to include that match the given patient id """
    include_list = []
    
    with open(keep_research_samples_file, newline='') as include:
        reader = csv.reader(include)
        for row in reader:
            if not row or not row[0].strip():
                continue
            sample_id = row[0].strip()
            if sample_id.startswith(f"{patient_id}-"):
                include_list.append(sample_id)
                print(f"Found sample in include list: {sample_id}")

    return include_list

def get_exclude_list(exclude_samples_file, patient_id):
    """ Read in the exclude csv file return samples to exclude that match the given patient id """
    exclude_list = []
    
    with open(exclude_samples_file, newline='') as exclude:
        reader = csv.reader(exclude)
        for row in reader:
            if not row or not row[0].strip():
                continue
            sample_id = row[0].strip()
            if sample_id.startswith(f"{patient_id}-"):
                exclude_list.append(sample_id)

    return exclude_list

def get_combined_patient_id(cmo_id, dmp_id):
    """ Create a combined patient id if both cmo and dmp id exist. Otherwise return whichever id exists """
    if cmo_id and dmp_id:
        return f"{cmo_id}_{dmp_id}"
    elif cmo_id:
        return cmo_id
    elif dmp_id:
        return dmp_id
    else:
        return ""

def save_to_json(sample_dict, patient_id):
    """ Save dictionary with patient and sample information to a JSON file """
    with open(f'{patient_id}_all_samples.json', "w") as out:
        json.dump(sample_dict[patient_id], out, indent=4)
    print(f'Saved {patient_id}_all_samples.json.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BAM paths.")
    parser.add_argument("--id_mapping_file", required=True)
    parser.add_argument("--keep_research_samples_file", required=False)
    parser.add_argument("--exclude_samples_file", required=False)
    parser.add_argument("--clinical_access_key_file", required=True)
    parser.add_argument("--clinical_impact_key_file", required=True)
    parser.add_argument("--research_access_bam_dir_template", required=True)
    parser.add_argument("--clinical_access_sample_regex_pattern", required=True)
    parser.add_argument("--clinical_impact_sample_regex_pattern", required=True)
    parser.add_argument("--research_access_manifest_file_template", required=True)
    parser.add_argument("--research_access_mutations_maf_template", required=True)
    args = parser.parse_args()

    get_all_samples(args.id_mapping_file, args.research_access_bam_dir_template, args.clinical_access_key_file, 
                   args.clinical_impact_key_file, args.keep_research_samples_file, args.exclude_samples_file, 
                   args.clinical_access_sample_regex_pattern, args.clinical_impact_sample_regex_pattern, 
                   args.research_access_manifest_file_template, args.research_access_mutations_maf_template)
