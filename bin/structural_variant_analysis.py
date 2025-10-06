import csv
import pandas as pd
import os
import argparse
import json
import re

##############################
# Helper functions
##############################

def parse_assay_info(sample_meta):
    """Extract source and assay from sample metadata."""
    assay_type = sample_meta.get("assay_type", "")
    parts = assay_type.split("_")
    source = parts[0] if len(parts) > 0 else ""
    assay = parts[1].upper() if len(parts) > 1 else ""
    return source, assay

SV_TYPE_MAP = {
    "INVERSION": "INV",
    "DELETION": "DEL",
    "INSERTION": "INS",
    "DUPLICATION": "DUP",
    "TRANSLOCATION": "TRA",
    # Add more if needed
}

def is_valid_path(path):
    if not os.path.isfile(os.path.realpath(path)):
        print(f"[WARNING] file not found: {path}.")
        return False
    return True

def load_patient_data(patient_json):
    with open(patient_json) as json_file:
        patient_data = json.load(json_file)
        cmo_id = patient_data['cmo_id']
        dmp_id = patient_data['dmp_id']
        combined_id = patient_data['combined_id']
    return patient_data, cmo_id, dmp_id, combined_id

def infer_research_access_sv_path(research_access_sv_template, cmo_id, sample_id):
    research_access_sv_path = research_access_sv_template.replace("{cmo_patient_id}", cmo_id).replace("{sample_id}", sample_id)
    if is_valid_path(research_access_sv_path):
        return research_access_sv_path
    else:
        return False

##############################
# Data extraction functions
##############################

def get_research_access_sv_calls(research_access_sv_template, patient_data, cmo_id, access_structural_variant_gene_list):
    research_access_sv_calls = []
    for sample_id, sample_data in patient_data["samples"].items():
        if sample_data.get('assay_type') == "research_access" and sample_data.get('tumor_normal') == "tumor":
            source, assay = parse_assay_info(sample_data)
            research_access_sv_path = infer_research_access_sv_path(research_access_sv_template, cmo_id, sample_id)
            if research_access_sv_path:
                with open(research_access_sv_path, 'r') as research_access_sv_file:
                    research_access_sv_data = csv.DictReader(research_access_sv_file, delimiter='\t')
                    for row in research_access_sv_data:
                        if row['Gene1'] in access_structural_variant_gene_list or row['Gene2'] in access_structural_variant_gene_list:
                            # Initialize gene1 and gene2 with original values
                            gene1 = row["Gene1"]
                            gene2 = row["Gene2"]
                            
                            # Check if Fusion column exists and parse it for gene names
                            if "Fusion" in row and row["Fusion"].strip() != "":
                                fusion_text = row["Fusion"].strip()
                                # Look for pattern like "Protein Fusion: in frame {EML4:ALK}"
                                fusion_match = re.search(r'\{([^:]+):([^}]+)\}', fusion_text)
                                if fusion_match:
                                    gene1 = fusion_match.group(1).strip()
                                    gene2 = fusion_match.group(2).strip()
                            
                            research_access_sv_calls.append({
                                "sample_id": row["TumorId"],
                                "sv_type": row["SV_Type"],
                                "gene1": gene1, 
                                "gene2": gene2,
                                "chr1": row["Chr1"],
                                "pos1": row["Pos1"],
                                "chr2": row["Chr2"], 
                                "pos2": row["Pos2"],
                                "split_read_count": row["SplitReadSupport"], 
                                "paired_read_count": row["PairEndReadSupport"], 
                                "total_read_count": row["TumorReadCount"], 
                                "info": row["Fusion"],
                                "source": source,
                                "assay": assay
                            })
    return research_access_sv_calls

def get_clinical_sv_calls(clinical_sv_file, patient_data, access_structural_variant_gene_list):
    clinical_sv_calls = []
    if not is_valid_path(clinical_sv_file):
        return clinical_sv_calls
    
    # Map from sample_id to sample_meta for clinical samples
    sample_meta_map = {s_id: meta for s_id, meta in patient_data["samples"].items()
                       if meta.get("assay_type", "").startswith("clinical")}
    
    with open(clinical_sv_file, 'r') as f:
        clinical_sv = csv.DictReader(f, delimiter="\t")
        for row in clinical_sv:
            sample_id = row.get("Sample_ID", "")
            if sample_id in sample_meta_map:
                meta = sample_meta_map[sample_id]
                if row.get('Site1_Hugo_Symbol', '') in access_structural_variant_gene_list or \
                   row.get('Site2_Hugo_Symbol', '') in access_structural_variant_gene_list:
                    source, assay = parse_assay_info(meta)
                    clinical_sv_calls.append({
                        "sample_id": sample_id,
                        "sv_type": SV_TYPE_MAP.get(row.get("Class",""), row.get("Class","")),
                        "gene1": row.get("Site1_Hugo_Symbol", ""),
                        "gene2": row.get("Site2_Hugo_Symbol", ""),
                        "chr1": row.get("Site1_Chromosome", ""),
                        "pos1": row.get("Site1_Position", ""),
                        "chr2": row.get("Site2_Chromosome", ""), 
                        "pos2": row.get("Site2_Position", ""),
                        "split_read_count": row.get("Tumor_Split_Read_Count", ""), 
                        "paired_read_count": row.get("Tumor_Paired_End_Read_Count", ""), 
                        "total_read_count": None, # always None, per your original
                        "info": row.get("Event_Info", ""),
                        "source": source,
                        "assay": assay
                    })
    return clinical_sv_calls

##############################
# Main processing function
##############################

def generate_sv_table(patient_json, research_access_sv_template, clinical_sv_file, access_structural_variant_gene_list):
    patient_data, cmo_id, dmp_id, combined_id = load_patient_data(patient_json)
    sv_columns = ["sample_id", "patient_id", "cmo_patient_id", "dmp_patient_id", "sv_type", "gene1", "gene2", "chr1", "pos1", "chr2", "pos2",
                  "split_read_count", "paired_read_count", "total_read_count", "info", "source", "assay"]

    # Collect calls
    research_access_sv_calls = get_research_access_sv_calls(research_access_sv_template, patient_data, cmo_id, access_structural_variant_gene_list)
    clinical_sv_calls = get_clinical_sv_calls(clinical_sv_file, patient_data, access_structural_variant_gene_list)

    research_access_sv_calls_df = pd.DataFrame(research_access_sv_calls, columns=sv_columns)
    clinical_sv_calls_df = pd.DataFrame(clinical_sv_calls, columns=sv_columns)

    # Add patient info columns
    research_access_sv_calls_df["patient_id"] = combined_id
    research_access_sv_calls_df["cmo_patient_id"] = cmo_id
    research_access_sv_calls_df["dmp_patient_id"] = dmp_id
    clinical_sv_calls_df["patient_id"] = combined_id
    clinical_sv_calls_df["cmo_patient_id"] = cmo_id
    clinical_sv_calls_df["dmp_patient_id"] = dmp_id

    # Combine frames
    all_sv_calls_df = pd.concat([research_access_sv_calls_df, clinical_sv_calls_df], ignore_index=True)

    # ========== Added columns, derived BEFORE logic ==========

    # variant = gene1__gene2:sv_type
    for df in (research_access_sv_calls_df, clinical_sv_calls_df, all_sv_calls_df):
        df["variant"] = df["gene1"].astype(str).str.strip() + "__" + \
                        df["gene2"].astype(str).str.strip() + ":" + \
                        df["sv_type"].astype(str).str.strip()

    # vaf = (split_read_count + paired_read_count)/total_read_count if all are numbers
    def calc_vaf(row):
        try:
            split = float(row["split_read_count"]) if pd.notna(row["split_read_count"]) and str(row["split_read_count"]).strip() != "" else 0
            paired = float(row["paired_read_count"]) if pd.notna(row["paired_read_count"]) and str(row["paired_read_count"]).strip() != "" else 0
            total = float(row["total_read_count"]) if pd.notna(row["total_read_count"]) and str(row["total_read_count"]).strip() != "" else None
            if total and total > 0:
                return (split + paired) / total
            else:
                return None
        except Exception:
            return None

    all_sv_calls_df["vaf"] = all_sv_calls_df.apply(calc_vaf, axis=1)

    # ========== Signed Out computed LAST ==========
    # Use clinical_sv_calls_df['variant'] for set, after it's created correctly
    clinical_variants = set(clinical_sv_calls_df["variant"].dropna().str.strip().tolist())

    def is_blank_genes(row):
        return str(row["gene1"]).strip() == "" and str(row["gene2"]).strip() == ""

    def variant_signed_out(row):
        if is_blank_genes(row):
            return None
        return "Signed Out" if str(row["variant"]).strip() in clinical_variants else None

    all_sv_calls_df["variant_clinical_status"] = all_sv_calls_df.apply(variant_signed_out, axis=1)
  
    # Create helper function for numerical chromosome sorting
    def chr_sort_key(chr_str):
        """Convert chromosome string to numerical value for proper sorting"""
        chr_str = str(chr_str).upper().replace('CHR', '')
        if chr_str == 'X':
            return 23
        elif chr_str == 'Y':
            return 24
        elif chr_str == 'MT' or chr_str == 'M':
            return 25
        else:
            try:
                return int(chr_str)
            except ValueError:
                return 26  # For any other non-standard chromosomes

    # Create temporary sorting columns for numerical chromosome sorting
    all_sv_calls_df['_chr_1_sort'] = all_sv_calls_df['chr1'].apply(chr_sort_key)
    all_sv_calls_df['_chr_2_sort'] = all_sv_calls_df['chr2'].apply(chr_sort_key)
    all_sv_calls_df['_pos_1_sort'] = pd.to_numeric(all_sv_calls_df['pos1'], errors='coerce')
    all_sv_calls_df['_pos_2_sort'] = pd.to_numeric(all_sv_calls_df['pos2'], errors='coerce')

    # Sort the DataFrame by numerical chromosomes, positions, and sample_id
    all_sv_calls_df = all_sv_calls_df.sort_values(
        by=['_chr_1_sort', '_pos_1_sort', '_chr_2_sort', '_pos_2_sort', 'sample_id'], 
        ignore_index=True
    )

    # Drop the temporary sorting columns
    all_sv_calls_df = all_sv_calls_df.drop(columns=['_chr_1_sort', '_chr_2_sort', '_pos_1_sort', '_pos_2_sort'])


    return all_sv_calls_df

##############################
# CLI
##############################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate SV table for a patient.")
    parser.add_argument("--patient_json", required=True, help="Path to patient JSON file.")
    parser.add_argument("--research_access_sv_template", required=True)
    parser.add_argument("--clinical_sv_file", required=True)
    parser.add_argument("--access_structural_variant_gene_list", required=True)
    parser.add_argument("--output", required=True, help="Output CSV file.")
    args = parser.parse_args()

    if os.path.isfile(args.access_structural_variant_gene_list):
        with open(args.access_structural_variant_gene_list) as f:
            access_structural_variant_gene_list = [x.strip() for x in f if x.strip()]
    else:
        access_structural_variant_gene_list = args.access_structural_variant_gene_list.split(",")

    df = generate_sv_table(args.patient_json, args.research_access_sv_template, args.clinical_sv_file,
                           access_structural_variant_gene_list)
    output_file = args.output
    df.to_csv(output_file, index=False)
    print(f'{output_file} has been created.')
