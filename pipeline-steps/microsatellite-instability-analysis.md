# Microsatellite Instability analysis

## Steps

* **Infer MSI file paths:**\
  Identify and locate paths to research and clinical ACCESS & IMPACT MSI output files for all patient samples.
* **Aggregate and annotate MSI calls:**
  * Extract MSI results for research ACCESS samples.
  * Extract MSI results for clinical ACCESS samples from specified `clinical_access_msi_file`.
  * Extract MSI results for clinical IMPACT samples from specified `clinical_impact_msi_file`.
  * Compile all extracted MSI calls into a single table with standardized columns.
* **Output:** `{patient_id}.msi.csv`

## Output file

Output file is organized in:

```
{outdir}/
└── final/
    └── {patient_id}/
        └── {patient_id}.msi.csv
```
