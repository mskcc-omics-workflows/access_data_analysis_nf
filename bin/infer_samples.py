import json
import pandas as pd
import argparse
import csv
import re
import os

"""
Script to read in the CMO/DMP sample IDs and retrieve all associated 
research access, clinical access, and clinical impact samples.

Samples specified in the include and exclude files will be filtered as well.

Output is one JSON file per patient, containing all samples relevant to the patient under a combined cmo/dmp id.
"""

def get_all_samples(id_mapping_file, research_access_bam_dir_template, clinical_access_key_file, clinical_impact_key_file, include_samples_file, exclude_samples_file, clinical_access_sample_regex_pattern, clinical_impact_sample_regex_pattern):
    """ Main logic function to get all samples from the id mapping file, split them by patient, get the relevant samples, and save to JSON. """

    # Extract the cmo ids, dmp ids, and combined ids from the input file.
    id_list = get_id_mapping(id_mapping_file)

    # Go through each patient one by one
    for patient_data in id_list:
        
        # Get all the ids of the patient
        combined_id = patient_data["combined_id"]
        cmo_id = patient_data["cmo_id"]
        dmp_id = patient_data["dmp_id"]
        sample_dict = { combined_id: patient_data }

        # If the patient has a cmo id, find any samples that need to be included/excluded, then get all research samples
        if cmo_id: 
            cmo_include_list = get_include_list(include_samples_file, cmo_id)
            cmo_exclude_list = get_exclude_list(exclude_samples_file, cmo_id)
            find_research_samples(cmo_include_list, cmo_exclude_list, research_access_bam_dir_template, cmo_id, combined_id, sample_dict)
        
        # If the patient has a dmp id, find any samples that need to be included/excluded, then get all clinical samples
        if dmp_id:
            dmp_exclude_list = get_exclude_list(exclude_samples_file, dmp_id)
            find_clinical_samples(clinical_impact_key_file, clinical_access_key_file, dmp_id, combined_id, dmp_exclude_list, sample_dict, clinical_access_sample_regex_pattern, clinical_impact_sample_regex_pattern)
 
        # Save the final dictionary to json
        save_to_json(sample_dict, combined_id)
    
    if not id_list:
        print("No samples found in input file.")


def find_research_samples(include_list, exclude_list, research_access_bam_dir_template, cmo_id, combined_id, sample_dict):

    # Get the root path of the access bam directory for the patient by removing the sample/current 
    research_access_bam_dir_root = research_access_bam_dir_template.split("/{sample_id}")[0].replace("{cmo_patient_id}", cmo_id)

    # Use the folders in the directory to get the research sample names
    research_sample_names = [ sample_name for sample_name in os.listdir(research_access_bam_dir_root)
        # Only keep the research samples that have a "current" folder and at least one existing bam file
        if is_valid_research_sample(sample_name, research_access_bam_dir_root)]

    # include and exclude samples based on input lists
    filtered_research_samples = filter_research_samples(research_sample_names, include_list, exclude_list)

    # add each sample to the sample dictionary
    for sample_id in filtered_research_samples:
        sample_dict[combined_id]["samples"][sample_id] = {
            "sample_id": sample_id,
            # infer whether the sample is tumor or normal
            "tumor_normal": infer_tumor_normal(sample_id),
            "assay_type": "research_access",
            "anon_id": "NA"
            }

def is_valid_research_sample(sample_name, research_access_bam_dir):
    """ Checks if an access research sample has a current directory with at least one bam file """

    # Check if current directory exists
    current_path = os.path.join(research_access_bam_dir, sample_name, "current")
    if not os.path.isdir(current_path):
        print(f'{sample_name} missing current directory.')
        return False

    # Check for bam file
    if not any(f.endswith(".bam") for f in os.listdir(current_path)):
        print(f'{sample_name} has no bam files in current directory.')
        return False

    return True

def filter_research_samples(sample_list, include_list, exclude_list):
    """ Go through the list of research samples and adds from the include list and removes from the exclude list """

    # Add samples from the include list to the sample list if they aren't already there
    for include_sample in include_list:
        if include_sample not in sample_list:
            print(f'Including sample {include_sample}.')
            sample_list.append(include_sample)
        else:
            print(f'{include_sample} already included.')
    
    # Create a final sample list add add samples to it if they are not on the exclude list
    final_sample_list = []
    for sample in sample_list:
        if sample in exclude_list:
            print(f'Excluded sample: {sample}.')
        else:
            final_sample_list.append(sample)

    return final_sample_list

def infer_tumor_normal(sample_id):
    """ Determine whether sample is tumor or normal, assuming sample format of C-XXXX-NXXXX or P-XXXX-NXXXX """

    # Use the '-' to split the patient id
    if sample_id.split("-")[2].startswith("N"):
        sample_type = "normal"
    else:
        sample_type = "tumor"
    
    return sample_type

def find_clinical_samples(clinical_impact_key_file, clinical_access_key_file, dmp_id, combined_id, exclude_list, sample_dict, clinical_access_sample_regex_pattern, clinical_impact_sample_regex_pattern):
    
    # Search the impact key file for matching samples by dmp_id
    clinical_access_samples = search_dmp_key_file(clinical_access_key_file, "clinical_access", clinical_access_sample_regex_pattern, combined_id, dmp_id, exclude_list, sample_dict)
    clinical_impact_samples = search_dmp_key_file(clinical_impact_key_file, "clinical_impact", clinical_impact_sample_regex_pattern, combined_id, dmp_id, exclude_list, sample_dict)

def search_dmp_key_file(dmp_key_path, assay_type, regex_pattern, combined_id, dmp_id, exclude_list, sample_dict):
    # Go through each line in the key file and search for the dmp id and sample pattern
    with open(dmp_key_path, 'r') as key_file:
        for line in key_file:
            if re.search(rf"{dmp_id}-{regex_pattern}", line):
                # Get the sample name and anon id from the line in the key file
                sample_id, anon_id, sample_type = parse_key_line(line)
                # Skip the sample if it is in the exclude list, otherwise add the sample to the sample dictionary
                if sample_id in exclude_list:
                    print(f"Excluded sample: {sample_id}.")
                    continue
                sample_dict[combined_id]["samples"][sample_id] = { "sample_id": sample_id, "tumor_normal": sample_type, "assay_type": assay_type, "anon_id": anon_id }


# Function to get the sample_id, sample_type, and anon_id from a line in an impact or access key file
def parse_key_line(line):

    # Split the line by columns. The first column in the key file is the sample id and the second is the anon id
    cols = line.split(sep=',')
    sample_id = cols[0]
    anon_id = cols[1]

    # Use the sample id to infer if sample is tumor or normal 
    sample_type = infer_tumor_normal(sample_id)

    # Remove anon id suffix from access samples
    if "-XS" in sample_id:
        anon_id = re.sub(r'-standard$', '', anon_id)

    return sample_id, anon_id, sample_type

def get_id_mapping(id_mapping_file):
    """ Read in the cmo and dmp ids, generate the combined patient id, and save each patient id in a list. """
    id_list = []

    # Go through each row of the input id csv and extract the ids
    with open(id_mapping_file, newline='') as ids:
        reader = csv.DictReader(ids)
        for row in reader:
            cmo_id = row['cmo_patient_id'].strip()
            dmp_id = row['dmp_patient_id'].strip()
            # Create a combined cmo/dmp id
            combined_id = get_combined_patient_id(cmo_id, dmp_id)
            # Add an entry with all the patient data to the list if the row is not empty
            if cmo_id or dmp_id:
                id_list.append({ "combined_id": combined_id, "cmo_id": cmo_id, "dmp_id": dmp_id, "samples": {} })

    return id_list

def get_include_list(include_samples_file, patient_id):
    """ Read in the include csv file and return samples to include that match the given patient id """
    include_list = []
    with open(include_samples_file, newline='') as include:
        reader = csv.reader(include)
        # Go through each sample
        for row in reader:
            # Skip empty lines
            if not row or not row[0].strip():
                continue
            # Get the sample id
            sample_id = row[0].strip()
            # Check if the sample id belongs to the patient
            if sample_id.startswith(f"{patient_id}-"):
                include_list.append(sample_id)

    return include_list

def get_exclude_list(exclude_samples_file, patient_id):
    """ Read in the exclude csv file return samples to exclude that match the given patient id """
    exclude_list = []
    with open(exclude_samples_file, newline='') as exclude:
        reader = csv.reader(exclude)
        for row in reader:
            # Skip empty lines
            if not row or not row[0].strip():
                continue
            # Get the sample id
            sample_id = row[0].strip()
            # Check if the sample id belongs to the patient
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
    parser.add_argument("--include_samples_file", required=False)
    parser.add_argument("--exclude_samples_file", required=False)
    parser.add_argument("--clinical_access_key_file", required=True)
    parser.add_argument("--clinical_impact_key_file", required=True)
    parser.add_argument("--research_access_bam_dir_template", required=True)
    parser.add_argument("--clinical_access_sample_regex_pattern", required=True)
    parser.add_argument("--clinical_impact_sample_regex_pattern", required=True)
    args = parser.parse_args()

    get_all_samples(args.id_mapping_file, args.research_access_bam_dir_template, args.clinical_access_key_file, args.clinical_impact_key_file, args.include_samples_file, args.exclude_samples_file, args.clinical_access_sample_regex_pattern, args.clinical_impact_sample_regex_pattern)