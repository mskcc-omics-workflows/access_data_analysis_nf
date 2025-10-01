# Copy Number analysis

## Steps

* **Infer CNA file paths:**\
  Determine paths to clinical and research ACCESS & IMPACT CNA files for all patient samples.
* **Aggregate and annotate variants:**
  * Combine all variant calls across samples.
  * Assign `variant_clinical_status = 'Signed Out'` to variants observed in clinical samples from the patient.
    * `variant` is defined as `<Hugo_Symbol>_<cna_type>`.
    * `cna_type` is `AMP` if fold change (`fc`) > 0; `DEL` if `fc` < 0.
    * For variants called in ACCESS samples, set `Signed Out` status only if the `<Hugo_Symbol>_<cna_type>` matches a variant present in at least one of the clinical samples.
* **Apply Filtering and Document Reasons:**
  * Create a `filter` column indicating reasons for exclusion or annotation.
    * **`pval_filter`:** Variant has a p-value above the threshold in research ACCESS sample.
    * **`fc_filter_signedout`:** Signed Out variant does not meet fold change threshold for signed out events (default thresholds: -1.2, 1.2) in research ACCESS sample.
    * **`fc_filter_denovo`:** Variant is in the designated CNA gene list but does not meet fold change threshold for de novo events (default thresholds: -1.5, 1.5) in research ACCESS sample.
    * **`denovo_not_in_genelist`:** Variant is neither signed out nor in the ACCESS CNA gene list.
    * **`gene_not_in_access`:** Clinical IMPACT variant where the gene was not observed in any ACCESS result files.
  * If multiple filtering reasons apply, they are combined in the filter column, separated by semicolons.
  * Variants passing all filters are marked as `PASS`.
  * Variants that are not PASS in any of the patient samples are excluded from all results files.
* **Output Results:**
  * All CNV calls: `{patient_id}.cnv.csv`
  * PASS-filtered CNV calls: `{patient_id}.cnv.pass-filtered.csv`

## Output files

Output files are organized as:

```
{outdir}/
└── intermediate/
│   └── {patient_id}/
│       └── {patient_id}.cnv.csv
└── final/
    └── {patient_id}/
        └── {patient_id}.cnv.pass-filtered.csv
```
