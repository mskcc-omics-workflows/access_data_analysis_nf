# Microsatellite Instability analysis

## Steps

* **Infer MSI File Paths:**\
  Identify and locate paths to research-access and clinical MSI score files for all relevant patient samples.
* **Aggregate and Annotate MSI Calls:**
  * Extract MSI scores from research ACCESS samples.
  * Extract MSI scores from clinical ACCESS samples from specified `clinical_access_msi_file`.
  * Extract MSI scores from clinical IMPACT samples from specified `clinical_impact_msi_file`.
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
