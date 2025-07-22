
# infer_samples.py

The script processes CMO and DMP patient IDs from an input CSV file, finds associated research and clinical samples, filters them based on optional include/exclude lists, and saves the results in one JSON file per patient.

## Inputs

Three input files:
- id_mapping_file
   - The input samples CSV file has two columns `cmo_patient_id` and `dmp_patient_id`. Even if you are only providing one type of id, both headers are required. The script can handle any combination of ids, it does not have to be standard for each sample.  
- include_samples_file
   - Optional CSV file with a list of sample ids to include. No header is necessary. The file can be empty, but needs to exist for the pipeline.
- exclude_samples_file
   - Optional CSV file with a list of sample ids to exclude. No header is necessary. The file can be empty, but needs to exist for the pipeline. 

## Assumptions

- The input ID samples csv file has two columns named exactly `cmo_patient_id` and `dmp_patient_id`. The include csv and exclude csv exist.
- Research samples:
   - Research samples are stored in directories defined by the `research_access_dir_template`, which includes a placeholder `{cmo_patient_id}`.
   - Each research sample directory contains a subdirectory named `"current"` that contains the relevant BAM files (`*.bam`).
- Clinical samples:
   - There are multiple bams for each clinical ACCESS sample. The script searches for the `standard` bam in the key file and uses that entry to get the relevant sample information. This is declared in the `access_sample_regex_pattern` in the config file. 
   - The IMPACT clinical samples also use a regex pattern, `impact_sample_regex_pattern` that assumes the sample name contains `-IM or -IH`. These regex patterns can be edited in the config file. 
   - The first column of the key file is the sample id and the second column of the key file is the anon id. The function `parse_key_file` will need to change if the key file structure changes.  
- Tumor/normal:
   Assumes the normal samples follow a format of `X-XXXX-NXXXX`. The function `infer_tumor_normal` will need to change if the samples do not follow this structure. 

## Output JSON Structure

Each JSON file is named `<combined_patient_id>_all_samples.json` and contains a single key: the combined patient ID (e.g., `"C-12345_P-67890"`). The value is a dictionary structured as follows:

```json
{
  "C-12345_P-67890": {
    "combined_id": "C-12345_P-67890",
    "cmo_id": "C-12345",
    "dmp_id": "P-67890",
    "samples": {
      "sample_id_1": {
        "sample_id": "sample_id_1",
        "tumor_normal": "tumor",
        "assay_type": "research_access",
        "anon_id": "NA"
      },
      "sample_id_2": {
        "sample_id": "sample_id_2",
        "tumor_normal": "normal",
        "assay_type": "clinical_access",
        "anon_id": "anon123"
      }
      // ... additional samples
    }
  }
}

```

---

## How the Script Works

1. **Read IDs:**  
   Reads the CSV ID file and generates a list of patient records. Each record includes CMO ID, DMP ID, combined patient ID, and an dictionary to store samples.

2. **Process Each Patient:**  
   For each patient record:
   - Retrieve research samples based on the `research_access_dir_template`.
   - Retrieve clinical samples from DMP key files using regex patterns.
   - Apply inclusion and exclusion filters to both sample sets.
   - Add the final samples to the patient's samples dictionary.

3. **Save Output:**  
   Saves the samples dictionary for each patient as a JSON file named `<combined_patient_id>_all_samples.json`.

## Getting Research Samples

- The script constructs the path to research samples by replacing `{cmo_patient_id}` in the `research_access_dir_template`.
- It lists all sample directories under this path.
- Each sample directory must contain a `current` folder with at least one `.bam` file to be considered valid.
- Samples from the include list are added if missing; samples in the exclude list are removed.
- Each valid sample is added to the patient’s dictionary with:
  - `sample_id`
  - `tumor_normal` status inferred from the sample ID (`-N` = normal, otherwise tumor)
  - `assay_type` set as `"research_access"`
  - `anon_id` set as `"NA"`

## Getting Clinical Samples

- Clinical samples are found by searching the DMP key files (`dmp_access_key_path` for clinical access, `dmp_impact_key_path` for clinical impact).
- The script uses provided regex patterns to match sample IDs relevant to the patient’s DMP ID.
- Samples listed in the exclude list are skipped.
- Each matched sample is added with:
  - `sample_id`
  - `tumor_normal` inferred from the sample ID
  - `assay_type` set as either `"clinical_access"` or `"clinical_impact"`
  - `anon_id` parsed and cleaned (removes `-standard` suffix if present for access samples)

---