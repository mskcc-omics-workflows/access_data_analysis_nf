import json
import os

def load_patient_data(patient_json):
    with open(patient_json) as json_file:
        patient_data = json.load(json_file)
        cmo_id = patient_data['cmo_id']
        dmp_id = patient_data['dmp_id']
        combined_id = patient_data['combined_id']
    return patient_data, cmo_id, dmp_id, combined_id

def save_to_csv(df, patient_id, var_tag):
    output_file = f'{patient_id}_{var_tag}.csv'
    df.to_csv(output_file, index=False)
    print(f'{output_file} has been created.')

# can be used for tsv, maf, txt
def save_to_tsv(df, patient_id, var_tag, suffix):
    output_file = f'{patient_id}_{var_tag}.{suffix}'
    df.to_csv(output_file, sep="\t", index=False)
    print(f'{output_file} has been created.')

def is_valid_path(path):
    if not os.path.isfile(os.path.realpath(path)):
        print(f"[WARNING] file not found: {path}.")
        return False
    return True