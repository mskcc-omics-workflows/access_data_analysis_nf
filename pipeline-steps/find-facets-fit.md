# Find FACETS fit

* Selects optimal FACETS fit for each IMPACT tumor sample (prioritizing reviewed/QC-passed).
* Falls back to default fit if needed.
* Outputs tab-separated mapping of samples to FACETS fits.
* Marks missing cases with `"MISSING"`.
*   Output file saved in:&#x20;

    ```
    {outdir}/
    └── intermediate/
        └── {patient_id}/
           └── {patient_id}_facets_fit.txt
    ```
