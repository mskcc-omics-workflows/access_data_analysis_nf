# Structural Variant analysis

## Steps

* **Infer SV file paths:**\
  Determine paths to clinical and research ACCESS & IMPACT SV files for all patient samples.
* **Aggregate and annotate SV calls:**
  * Combine SV calls from research and clinical samples.
  * Map SV types for clinical samples using a standardized nomenclature (e.g., `INV`, `DEL`, `BND`, etc.).
  * Build a unique `variant` identifier such as `<gene_1>__<gene_2>:<sv_type>`.
  * For each variant observed in clinical samples, assign `variant_clinical_status = 'Signed Out'`.
  * For research ACCESS-only variants, keep only those variants that involves at least one gene in the specified ACCESS SV gene list.
* Output file: `{patient_id}.sv.csv`

## Output file

Output file is organized in:

```
{outdir}/
└── final/
    └── {patient_id}/
        └── {patient_id}.sv.csv
```
