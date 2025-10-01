# Biometrics

* Runs steps from the python package [biometrics](https://cmo-ci.gitbook.io/biometrics/).
* Create biometrics input file using the patient json and standard BAM path templates
* _**Extraction**_ to get pileup and coverage information from the BAM files
* _**Genotype**_ for finding sample matches and mismatches
* _**Sex mismatch**_ determination
* Based on the genotype and sex mismatch output files to create a CSV _**summary**_ file.
  * genotype\_qc\_status is marked FAIL if there are any unexpected mismatches, and sex\_mismatch\_qc\_status is marked FAIL if any of the samples has a sex mismatch.



Output files are organized as:

```
{outdir}/
└── intermediate/
│   └── biometrics/
│      └── {patient_id}/
│           ├── {patient_id}.biometrics_input.csv
│           ├── {patient_id}.genotype_comparison.csv
│           ├── {patient_id}.sex_mismatch.csv
│           └── extract_db/
└── final/
    └── {patient_id}/
        └── {patient_id}.biometrics_summary.csv
```
