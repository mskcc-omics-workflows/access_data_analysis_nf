# Structural Variant analysis

## Steps

* **Infer SV File Paths:**\
  Determine paths to clinical and research-access SV files for all patient samples.
* **Aggregate and Annotate SV Calls:**
  * Combine SV calls from research and clinical sources across all samples.
  * Map SV types for clinical samples using a standardized nomenclature (e.g., `INV`, `DEL`, `BND`, etc.).
  * Build a unique `variant` identifier such as `<gene_1>__<gene_2>:<sv_type>`.
  * For each variant observed in clinical samples, assign `variant_clinical_status = 'Signed Out'`.
  * For research ACCESS calls, annotate only those variants involving genes in the specified ACCESS SV gene list.
* Output file: `{patient_id}.sv.csv`

## Output file

Output file is organized in:

```
{outdir}/
└── final/
    └── {patient_id}/
        └── {patient_id}.sv.csv
```
