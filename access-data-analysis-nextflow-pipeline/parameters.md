# Parameters

All parameters are set in `nextflow.config` file, and their schema specified in `nextflow_schema.json`.  Juno-specific paths specified in `conf/juno.config`.&#x20;

TO DO: Create a config file for iris: `conf/iris.config`.&#x20;

## Input and Output parameters

* **input** (_required_): CSV file with columns: `cmo_patient_id`, `dmp_patient_id`, `sex` of patients to include.
  * At least one of cmo\_patient\_id or dmp\_patient\_id is required.
  * sex is required, and must be either F or M.
* **keep\_research\_samples\_file** (_optional_): File containing _**CMO sample IDs**_ to include (no header, one sample id per row).
  * If provided, CMO samples (both tumor and normal) not in the list will be ignored.
  * _Use this to provide the same sample list you would specify in the R pipeline masterfile._
  * Use `null` to not restrict CMO samples to a list.
* **exclude\_samples\_file** (_optional_): File containing _**CMO sample IDs**_ and _**DMP sample IDs**_ to exclude (no header, one sample id per row).
  * Exclusions override inclusions for CMO sample IDs.
  * Use `null` to not exclude any samples.
* **outdir** (_required_): Output directory path.

&#x20;
