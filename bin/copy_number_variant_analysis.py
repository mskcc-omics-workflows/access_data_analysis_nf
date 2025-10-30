import csv
import pandas as pd
import os
import argparse
import json
import numpy as np

def load_patient_data(patient_json):
    """Load patient metadata and sample info from a JSON file."""
    with open(patient_json) as json_file:
        patient_data = json.load(json_file)
    return (
        patient_data,
        patient_data['cmo_id'],
        patient_data['dmp_id'],
        patient_data['combined_id'],
    )

def parse_assay_info(sample_meta):
    """Extract source and assay from sample metadata."""
    assay_type = sample_meta.get("assay_type", "")
    parts = assay_type.split("_")
    source = parts[0] if len(parts) > 0 else ""
    assay = parts[1].upper() if len(parts) > 1 else ""
    return source, assay

def get_clinical_signed_out_cnas(patient_data, dmp_id, clinical_cna_file, combined_id, access_gene_list_v1, access_gene_list_v2):
    clinical_calls = []
    clinical_gene_event_set = set()
    if dmp_id:
        with open(clinical_cna_file, 'r') as f:
            reader = csv.reader(f, delimiter="\t")
            data = list(reader)
            header = data[0]
            for col_index, sample_id in enumerate(header[1:], start=1):
                if str(dmp_id) in sample_id:
                    for row in data[1:]:
                        gene = row[0]
                        cna_val = row[col_index]
                        if cna_val and cna_val != '0':
                            cna_val = float(cna_val)
                            cna_type = "AMP" if cna_val > 0 else "DEL" if cna_val < 0 else ""
                            if not cna_type:
                                continue
                            clinical_gene_event_set.add((gene, cna_type))
                            sample_meta = patient_data["samples"].get(sample_id, {})
                            source, assay = parse_assay_info(sample_meta)
                            clinical_calls.append({
                                "Hugo_Symbol": gene,
                                "cna_type": cna_type,
                                "sample_id": sample_id,
                                "patient_id": combined_id,
                                "cmo_patient_id": "",
                                "dmp_patient_id": dmp_id,
                                "fold_change": np.nan,
                                "p_val": "",
                                "filter": "PASS",
                                "assay": assay,
                                "source": source
                            })
    return pd.DataFrame(clinical_calls), clinical_gene_event_set

def infer_research_cna_path(template, cmo_id, sample_id):
    path = template.replace("{cmo_patient_id}", cmo_id).replace("{sample_id}", sample_id)
    if not os.path.isfile(os.path.realpath(path)):
        print(f"[WARNING] CNA file not found: {path}")
        return None
    return path

def process_research_access_calls(
    patient_data, cmo_id, combined_id, access_gene_list_v1, access_gene_list_v2,
    clinical_gene_event_set, research_access_cna_template,
    pval_threshold, fc_denovo_amp, fc_denovo_del, fc_signedout_amp, fc_signedout_del
):
    research_calls = []
    for sample_id, sample_data in patient_data["samples"].items():
        access_version = sample_data.get("access_version")
        if sample_data.get("assay_type") != "research_access" or sample_data.get("tumor_normal") != "tumor":
            continue
        cna_path = infer_research_cna_path(research_access_cna_template, cmo_id, sample_id)
        if not cna_path:
            continue
        source, assay = parse_assay_info(sample_data)
        with open(cna_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene = row.get("region", "")
                try:
                    fc = float(row['fc'])
                    p_val = float(row["p.adj"])
                except Exception:
                    fc = np.nan
                    p_val = np.nan
                cna_type = "AMP" if fc > 0 else "DEL" if fc < 0 else ""
                filter_reasons = []
                gene_lists = []
                if not np.isnan(p_val) and p_val > pval_threshold:
                    filter_reasons.append("pval_filter")
                is_signedout_event = (gene, cna_type) in clinical_gene_event_set
                if cna_type in ("AMP", "DEL"):
                    if is_signedout_event:
                        if (cna_type == "AMP" and fc < fc_signedout_amp) or \
                           (cna_type == "DEL" and fc > fc_signedout_del):
                            filter_reasons.append("fc_filter_signedout")
                    else:
                        if access_version == "XS1" and gene in access_gene_list_v1:
                            if (cna_type == "AMP" and fc < fc_denovo_amp) or \
                               (cna_type == "DEL" and fc > fc_denovo_del):
                                filter_reasons.append("fc_filter_denovo")
                        if access_version == "XS2" and gene in access_gene_list_v2:
                            if (cna_type == "AMP" and fc < fc_denovo_amp) or \
                               (cna_type == "DEL" and fc > fc_denovo_del):
                                filter_reasons.append("fc_filter_denovo")
                        else:
                            filter_reasons.append("denovo_not_in_genelist")
                else:
                    filter_reasons.append("unknown_event")
                filter_str = ";".join(filter_reasons) if filter_reasons else "PASS"
                research_calls.append({
                    "Hugo_Symbol": gene,
                    "cna_type": cna_type,
                    "sample_id": sample_id,
                    "patient_id": combined_id,
                    "cmo_patient_id": cmo_id,
                    "dmp_patient_id": "",
                    "fold_change": fc,
                    "p_val": p_val,
                    "filter": filter_str,
                    "assay": assay,
                    "source": source
                })
    return pd.DataFrame(research_calls)

def final_filter(df):
    """
    Mark 'gene_not_in_access' in filter if Hugo_Symbol not present in ACCESS samples.
    Do not drop any rows for this reason.
    """
    # Only keep genes with at least one sample where filter == PASS
    keep_genes = df.groupby('Hugo_Symbol')['filter'].apply(lambda x: (x == 'PASS').any())
    keep_genes = keep_genes[keep_genes].index
    df = df[df['Hugo_Symbol'].isin(keep_genes)].copy()

    access_genes = set(df[df['assay'] == "ACCESS"]['Hugo_Symbol'])
    mask = ~df['Hugo_Symbol'].isin(access_genes)

    # Rows where filter is empty or PASS
    pass_idx = mask & df['filter'].isin(['', 'PASS'])
    df.loc[pass_idx, 'filter'] = 'gene_not_in_access'

    # Rows where filter already has other entries
    not_pass_idx = mask & ~df['filter'].isin(['', 'PASS'])
    # Only append if 'gene_not_in_access' is not already present
    df.loc[not_pass_idx, 'filter'] = df.loc[not_pass_idx, 'filter'].apply(
        lambda x: x if 'gene_not_in_access' in x else x + ';gene_not_in_access'
    )

    df = df.sort_values(by=["Hugo_Symbol", "cna_type", "sample_id"]).reset_index(drop=True)
    return df

def final_filter_old(df):
    """
    Mark 'gene_not_in_access' in filter if Hugo_Symbol not present in ACCESS samples.
    Do not drop any rows for this reason.
    """
    # Only keep genes with at least one sample where filter == PASS
    keep_genes = df.groupby('Hugo_Symbol')['filter'].apply(lambda x: (x == 'PASS').any())
    keep_genes = keep_genes[keep_genes].index
    df = df[df['Hugo_Symbol'].isin(keep_genes)].copy()

    access_genes = set(df[df['assay'] == "ACCESS"]['Hugo_Symbol'])
    mask = ~df['Hugo_Symbol'].isin(access_genes)
    # Update filter entries for genes not in access: append or set
    pass_idx = mask & df['filter'].isin(['', 'PASS'])
    df.loc[pass_idx, 'filter'] = 'gene_not_in_access'
    not_pass_idx = mask & ~df['filter'].isin(['', 'PASS'])
    df.loc[not_pass_idx, 'filter'] = df.loc[not_pass_idx, 'filter'] + ';gene_not_in_access'
    df = df.sort_values(by=["Hugo_Symbol", "cna_type", "sample_id"]).reset_index(drop=True)
    return df

def add_variant_columns(df, clinical_gene_event_set):
    """Add variant and variant_clinical_status columns."""
    df["variant"] = df["Hugo_Symbol"] + "_" + df["cna_type"]
    df["variant_clinical_status"] = df.apply(
        lambda row: "Signed Out" if (row["Hugo_Symbol"], row["cna_type"]) in clinical_gene_event_set else "",
        axis=1
    )
    return df

def save_to_csv(df, output_file):
    df.to_csv(output_file, index=False)
    print(f"{output_file} has been created.")

def main(args):
    access_gene_list_v1 = args.access_copy_number_gene_list_v1.split(",")
    access_gene_list_v2 = args.access_copy_number_gene_list_v2.split(",")
    patient_data, cmo_id, dmp_id, combined_id = load_patient_data(args.patient_json)
    clinical_cnas, clinical_gene_event_set = get_clinical_signed_out_cnas(
        patient_data, dmp_id, args.clinical_cna_file, combined_id, access_gene_list_v1, access_gene_list_v2
    )
    research_cnas = process_research_access_calls(
        patient_data, cmo_id, combined_id, access_gene_list_v1, access_gene_list_v2,
        clinical_gene_event_set, args.research_access_cna_template,
        args.p_value_threshold, args.fc_denovo_amp, args.fc_denovo_del,
        args.fc_signedout_amp, args.fc_signedout_del
    )
    all_calls = pd.concat([clinical_cnas, research_cnas], ignore_index=True)
    all_calls = final_filter(all_calls)
    all_calls = add_variant_columns(all_calls, clinical_gene_event_set)
    all_calls = all_calls[
        ["Hugo_Symbol", "cna_type", "variant", "sample_id", "patient_id", "cmo_patient_id",
         "dmp_patient_id", "fold_change", "p_val", "variant_clinical_status", "filter", "assay", "source"]
    ]
    save_to_csv(all_calls, args.output)

    filtered_df = all_calls[all_calls["filter"].isin(["", "PASS"])]
    
    if not args.output_final:
        root, ext = os.path.splitext(args.output)
        args.output_final = f"{root}.pass-filtered.csv"
    save_to_csv(filtered_df, args.output_final)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compile clinical and research CNA calls for a patient.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--patient_json", required=True,
                        help="Path to patient JSON file containing metadata and sample info.")
    parser.add_argument("--research_access_cna_template", required=True,
                        help='Template for research access CNA file path. Use {cmo_patient_id} and {sample_id} as placeholders.')
    parser.add_argument("--clinical_cna_file", required=True,
                        help="TSV file with clinical signed-out CNA calls.")
    parser.add_argument("--access_copy_number_gene_list_v1", required=True,
                        help="Comma-separated list of genes eligible for ACCESS CNA reporting.")
    parser.add_argument("--access_copy_number_gene_list_v2", required=True,
                        help="Comma-separated list of genes eligible for ACCESS CNA reporting.")
    parser.add_argument("--p_value_threshold", type=float, default=0.05,
                        help="P-value threshold for filtering research ACCESS CNAs.")
    parser.add_argument("--fc_denovo_amp", type=float, default=1.5,
                        help="Fold-change threshold for de-novo amplifications (not previously signed out).")
    parser.add_argument("--fc_denovo_del", type=float, default=-1.5,
                        help="Fold-change threshold for de-novo deletions (not previously signed out).")
    parser.add_argument("--fc_signedout_amp", type=float, default=1.2,
                        help="Fold-change threshold for amplifications already signed out clinically.")
    parser.add_argument("--fc_signedout_del", type=float, default=-1.2,
                        help="Fold-change threshold for deletions already signed out clinically.")
    parser.add_argument("--output", required=True, help="Path to CSV file with all CNA calls (both PASS and non-PASS).")
    parser.add_argument("--output_final", required=False, help="Path to CSV file with final filtered CNA calls (PASS only). If not provided, script will write to <output_prefix>.pass-filtered.csv")
    args = parser.parse_args()
    main(args)
