# Setup and Usage

## Dependencies

* Nextflow ≥ 24.04.2
* Python ≥ 3.9
* pandas
* numpy

## Example Command

```
nextflow run main.nf -c nextflow.config -profile conda,juno_accessv1
```

## Parameters

All parameters are set in `nextflow.config` file, and their schema specified in `nextflow_schema.json`. &#x20;

## Input and Output parameters

* **`input`** (_required_): CSV file with columns: `cmo_patient_id`, `dmp_patient_id`, `sex` of patients to include.
  * At least one of cmo\_patient\_id or dmp\_patient\_id is required.
  * sex is required, and must be either F or M.
* **`keep_research_samples_file`** (_optional_): File containing _**CMO sample IDs**_ to include (no header, one sample id per row).
  * If provided, CMO samples (both tumor and normal) not in the list will be ignored.
  * _Use this to provide the same sample list you would specify in the R pipeline masterfile._
  * Use `null` to not restrict CMO samples to a list.
* **`exclude_samples_file`** (_optional_): File containing _**CMO sample IDs**_ and _**DMP sample IDs**_ to exclude (no header, one sample id per row).
  * Exclusions override inclusions for CMO sample IDs.
  * Use `null` to not exclude any samples.
* **`outdir`** (_required_): Output directory path.

## Profiles

* **`conda`** (currently required) – Creates and uses conda environments for the different pipeline steps.
* **`juno_accessv1`** (currently default)
  * Uses file paths specific to terra as defined in conf/juno\_accessv1.config
  * Uses DONOR22-TP in the small variants MAF file path
  * Uses CNV gene list specific to ACCESS v1
* **`juno_accessv2`**&#x20;
  * Uses file paths specific to terra as defined in conf/juno\_accessv2.config
  * Uses Donor19F21c2206-TP01  in the small variants MAF file path
  * Uses CNV gene list specific to ACCESS v2.
  * _\*Note: For biometrics, it uses the ACCESS v1 tiling SNPs file, as the larger number of tiling SNPs in ACCESS v2 significantly makes biometrics extraction step more computationally intensive and time-consuming._
* TO DO: Create `iris_accessv1` and `iris_accessv2` config files and profiles.&#x20;
