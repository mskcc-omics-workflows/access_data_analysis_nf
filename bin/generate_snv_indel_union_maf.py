import csv
import pandas as pd
import os
import argparse
import json
from collections import defaultdict
import glob

# Define the columns to be included in the MAF files
MAF_COLUMNS = [
    'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
    'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
    'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode',
    't_ref_count', 't_alt_count', 'n_ref_count', 'n_alt_count',
    'Variant_Classification', 'HGVSp', 'HGVSp_Short', 'HGVSc'
]

# Additional columns we'll add to the output
ADDITIONAL_COLUMNS = ['Clinical', 'Called_In']

# Variant identity columns - annotations may differ so exclude them from the key
VARIANT_KEY_COLS = [
    'Chromosome', 'Start_Position', 'End_Position',
    'Reference_Allele', 'Tumor_Seq_Allele2'
]

def get_all_calls(patient_json, research_access_mutations_maf_template, dmp_mutations_file):
    """
    Load patient data, get all mutation calls (research and clinical), merge and filter them,
    track their sources, then write results to a single MAF file.
    """
    # Load in the patient JSON
    patient_data = load_patient_data(patient_json)
    combined_id = patient_data.get('combined_id')
    # Get the research and clinical calls from corresponding MAF files
    research_calls, research_sample_map, research_assay_types = get_research_access_mutations(patient_data, research_access_mutations_maf_template)
    clinical_calls, clinical_sample_map, clinical_assay_types = get_clinical_mutations(patient_data, dmp_mutations_file)
    # Merge and filter the calls based on the exclude gene and classification lists
    all_small_calls, mutation_to_samples, sample_assay_types = merge_calls_with_tracking(
        research_calls, clinical_calls,
        research_sample_map, clinical_sample_map,
        research_assay_types, clinical_assay_types
    )

    # Write the final filtered output file with Called_In column
    final = create_maf(all_small_calls, mutation_to_samples, sample_assay_types, combined_id)
    return final

def parse_mutation_file(mutations_file, assay_type, dmp_id=None):
    """
    Parse maf files and extract all MAF columns needed.
    Returns:
    - list of mutations
    - dictionary mapping mutation identifiers to sample identifiers
    - dictionary mapping samples to their assay types
    """
    mutations = []
    sample_tracking = defaultdict(list)
    sample_assay_types = {}

    # Return empty results if the maf file doesn't exist
    if not os.path.exists(mutations_file):
        print(f"File does not exist: {mutations_file}")
        return mutations, sample_tracking, sample_assay_types

    try:
        with open(mutations_file, 'r') as maf:
            # For clinical mafs, exclude the metadata lines
            if assay_type == "clinical":
                maf_data = (line for line in maf if "sequenced_samples:" not in line)
            else:
                maf_data = maf

            reader = csv.DictReader(maf_data, delimiter='\t')
            for row in reader:
                # Exclude germline mutations from both assays
                if row.get('Mutation_Status') == 'GERMLINE':
                    continue
                # For research, remove variants that did not pass QC
                if assay_type == "research" and not row['Status'] == "":
                    continue
                # For clinical, only look at lines for the relevant dmp_id
                if assay_type == "clinical" and dmp_id and str(dmp_id) not in row['Tumor_Sample_Barcode']:
                    continue

                try:
                    # Create a variant dict with all possible columns
                    variant = {}
                    for col in MAF_COLUMNS:
                        variant[col] = row.get(col, '')

                    # Convert numeric fields
                    for numeric_field in ['Start_Position', 'End_Position']:
                        if variant[numeric_field]:
                            variant[numeric_field] = int(variant[numeric_field])

                    # Create a unique mutation identifier and track which sample it came from
                    mutation_key = (
                        f"{variant['Chromosome']}_{variant['Start_Position']}_"
                        f"{variant['End_Position']}_{variant['Reference_Allele']}_"
                        f"{variant['Tumor_Seq_Allele2']}"
                    )
                    sample_id = row['Tumor_Sample_Barcode']
                    
                    # Assign source priority for keeping annotations when called by multiple samples:
                    # 0 = clinical (highest priority), 1 = research (lowest priority)
                    variant['_source_priority'] = 0 if assay_type == "clinical" else 1
                    variant['_sample_id'] = sample_id

                    mutations.append(variant)
                    sample_tracking[mutation_key].append(sample_id)
                    # Record the assay type for this sample
                    sample_assay_types[sample_id] = assay_type

                except Exception as e:
                    print(f"Error parsing row in {mutations_file}: {e}")

    except Exception as e:
        print(f"Error opening or processing file {mutations_file}: {e}")

    return mutations, sample_tracking, sample_assay_types

def get_research_access_mutations(patient_data, research_access_mutations_maf_template):
    """
    Find every research sample in the patient_data, and parse the maf file for each sample.
    """
    cmo_id = patient_data.get('cmo_id')
    research_mutations = []
    combined_sample_tracking = defaultdict(list)
    combined_assay_types = {}

    if not cmo_id:
        return research_mutations, combined_sample_tracking, combined_assay_types

    # Find each research sample in the patient_data
    for sample_id, sample_data in patient_data.get("samples", {}).items():
        if sample_data.get('assay_type') == "research_access" and sample_data.get('tumor_normal') == "tumor":
            access_version = sample_data.get('access_version')
            donor_id = sample_data.get('donor_id')
            maf_path = (
                research_access_mutations_maf_template
                .replace("{cmo_patient_id}", cmo_id)
                .replace("{sample_id}", sample_id)
                .replace("{donor_id}", donor_id)
            )
            # Add the variants from the maf file
            mutations, sample_tracking, sample_assay_types = parse_mutation_file(maf_path, "research")
            research_mutations.extend(mutations)
            # Merge tracking information
            for mut_key, samples in sample_tracking.items():
                combined_sample_tracking[mut_key].extend(samples)
            combined_assay_types.update(sample_assay_types)

    return research_mutations, combined_sample_tracking, combined_assay_types

def get_clinical_mutations(patient_data, dmp_mutations_file):
    """ Get clinical mutations from the given dmp file. """
    dmp_id = patient_data.get('dmp_id')
    if not dmp_id or not os.path.exists(dmp_mutations_file):
        return [], {}, {}
    return parse_mutation_file(dmp_mutations_file, "clinical", dmp_id)

def merge_calls_with_tracking(research_calls, clinical_calls, research_sample_map,
                               clinical_sample_map, research_assay_types, clinical_assay_types):
    """
    Combine research and clinical calls and build a mapping of mutations to samples.
    Clinical annotations take priority over research annotations for the same variant.
    Within each source type, the first sample encountered (by name, descending) is used.
    """
    all_calls = research_calls + clinical_calls

    if not all_calls:
        all_calls_df = pd.DataFrame(columns=MAF_COLUMNS)
        return all_calls_df, {}, {}

    all_calls_df = pd.DataFrame(all_calls)

    # Build combined mutation-to-samples mapping
    mutation_to_samples = defaultdict(list)
    for mut_key, samples in research_sample_map.items():
        mutation_to_samples[mut_key].extend(samples)
    for mut_key, samples in clinical_sample_map.items():
        mutation_to_samples[mut_key].extend(samples)

    # Merge the assay type dictionaries
    sample_assay_types = {**research_assay_types, **clinical_assay_types}

    # Sort so that:
    # 1. Clinical rows come before research rows (_source_priority ascending: 0 before 1)
    # 2. Within each source type, sort by sample name descending as a tiebreaker
    all_calls_df = all_calls_df.sort_values(
        by=['_source_priority', '_sample_id'],
        ascending=[True, False]
    )

    # Drop duplicates keeping the first row per variant identity, which is now the
    # highest-priority annotation after sorting above
    all_calls_df = all_calls_df.drop_duplicates(subset=VARIANT_KEY_COLS, keep='first')

    # Drop internal helper columns before returning
    all_calls_df = all_calls_df.drop(columns=['_source_priority', '_sample_id'], errors='ignore')

    return all_calls_df, mutation_to_samples, sample_assay_types

def is_clinical_mutation(mutation_key, mutation_to_samples, sample_assay_types):
    """Check if a mutation appears in any clinical sample"""
    for sample in mutation_to_samples.get(mutation_key, []):
        if sample_assay_types.get(sample) == "clinical":
            return True
    return False

def create_maf(calls_df, mutation_to_samples, sample_assay_types, patient_id):
    """
    Save call information to a tab delimited file with Clinical status and Called_In columns,
    preserving all original columns.
    """
    if calls_df.empty:
        print(f"Warning: No mutations to write for {patient_id}")
        calls_df = pd.DataFrame(columns=MAF_COLUMNS + ADDITIONAL_COLUMNS)
    else:
        # Generate mutation keys for lookup
        mutation_keys = calls_df.apply(
            lambda row: (
                f"{row['Chromosome']}_{row['Start_Position']}_{row['End_Position']}_"
                f"{row['Reference_Allele']}_{row['Tumor_Seq_Allele2']}"
            ),
            axis=1
        )
        # Create Clinical column
        calls_df['Clinical'] = mutation_keys.apply(
            lambda key: "Signed Out" if is_clinical_mutation(key, mutation_to_samples, sample_assay_types) else ""
        )
        # Create Called_In column with semicolon-separated list of samples
        calls_df['Called_In'] = mutation_keys.apply(
            lambda key: ";".join(sorted(set(mutation_to_samples.get(key, []))))
        )

    return calls_df

def load_patient_data(patient_json):
    try:
        with open(patient_json) as json_file:
            return json.load(json_file)
    except Exception as e:
        print(f"Error loading patient data from {patient_json}: {e}")
        return {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate all called mutations MAF with sample tracking.")
    parser.add_argument("--patient_json", required=True, help="Path to patient JSON file")
    parser.add_argument("--research_access_mutations_maf_template", required=True,
                        help="Template path to research MAF files")
    parser.add_argument("--dmp_mutations_file", required=True, help="Path to clinical mutations file")
    parser.add_argument("--output", required=True, help="Output file")
    args = parser.parse_args()

    df = get_all_calls(
        args.patient_json,
        args.research_access_mutations_maf_template,
        args.dmp_mutations_file
    )
    output_file = args.output
    df.to_csv(output_file, index=False, sep="\t")
    print(f'{output_file} has been created with all columns preserved.')
