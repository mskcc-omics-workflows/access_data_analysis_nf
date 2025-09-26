# Pipeline Steps

## Research and Clinical ACCESS and IMPACT sample inference

* Parses ACCESS research BAM directories and key files for ACCESS/IMPACT clinical samples.
* Uses keep and exclude lists provided by user to keep or exclude samples.
* Sets combined patient ID as `{cmo_patient_id}_{dmp_patient_id}`
  * \[or `cmo_patient_id` or `dmp_patient_id` if other ID not present]
* Output files: `{patient_id}.json`.

The subsequent steps are run separately for each patient {patient\_id}.
