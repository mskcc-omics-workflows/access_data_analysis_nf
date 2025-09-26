# SNV/INDEL analysis

## 1. Union MAF generation

* Collects variants from ACCESS & IMPACT samples.
* Produces unified MAF file with all variants.
* _**\[NEW]**_ Does not filter any variants at this point
* _**\[NEW]**_ Adds columns:
  * **Called\_In**: Samples containing the variant (semicolon-separated).
  * **Clinical**: Marks variants from clinical samples as `"Signed Out"`
* _Output file: `{patient_id}_all_small_calls.maf`_

## 2. Input generation for genotype variants

* Builds genotyping input with inferred BAM paths:
  * Duplex & Simplex for ACCESS tumors.
  * Unfiltered for ACCESS normals.
  * Standard for IMPACT.
* Uses union MAF from previous step
* _Output file: `{patient_id}_genotyping_input.tsv`_

## 3. **Genotype Variants**

* Runs `genotype_variants small_variants multiple-samples`.
* Deduplication enabled.
* Requires `GetBaseCountsMultiSample` (path specified in the juno config file).
* Output files: `{patient_id}/genotyped_mafs/*.maf`

## 4. Aggregate variant allele counts across samples

* Produces per-sample allele counts + VAF for each variant.
* If there are M variants and N samples, there will be MxN rows in this table.
* Uses fragment counts for alt and total counts for ACCESS tumor samples and raw counts for ACCESS normal and IMPACT samples.
* Fragment counts are from the following BAM files:
  * ACCESS tumors: SIMPLEX + DUPLEX
  * ACCESS normals: UNFILTERED
  * IMPACT: STANDARD
* _**\[NEW]**_ Assigns call status per sample:
  * "called": variant was called in that sample
  * "low\_coverage": total count is below coverage threshold (default - ACCESS: 100x, IMPACT: 50x)
  * "genotyped": previous conditions not met and VAF > 0
  * "": (empty) none of the above conditions met
* Output file: `{patient_id}-SNV-INDEL.allele_counts.csv`

## 5. Annotate with Hotspot and CH lists

* Annotates against hotspot and clonal hematopoiesis lists (specified in config file).
* Adds **Hotspot** and **CH** columns ("yes" if matched, empty if not).
* Output file: `{patient_id}-SNV-INDEL.allele_counts.hotspot_ch.csv`

## 6. **Filter variants&#x20;**_**\[MODIFIED]**_

* _**\[NEW]**_ Adds 'filter' column to variant table, marking variants for exclusion based on:
  * Presence in excluded genes list ('excluded\_gene')
  * Matching excluded variant classifications ('excluded\_classification')
  * _**\[NEW]**_ Low coverage (low\_coverage) in ALL (tumor and normal, research and clinical) ACCESS samples ('low\_access\_cov')
  * _**\[MODIFIED]**_ For non-signed out variants if none of the ACCESS tumor samples meet minimum duplex alt threshold (default threshold: 3 for hotspot mutation, 5 for non-hotspot mutations) ('max\_duplex\_alt\_count')
  * _**\[MODIFIED]**_ For non-signed out variants if the ratio of max VAF of tumor samples to max VAF of normal samples if less than threshold (default: 2) ('low\_tumor\_to\_normal\_vaf\_ratio')
  * Multiple filter reasons are combined with semicolons (e.g., "excluded\_gene;low\_access\_cov")
* Creates three output files:
  1. _**\[NEW]**_ Full variant list with annotations. `{patient_id}-SNV-INDEL.allele_counts.hotspot_ch.filter.csv`
  2. Filtered PASS-only list. `{patient_id}.snv_indel.pass-filtered.csv`
  3. A table format with one row per variant of the pass-filtered file. `{patient_id}.snv_indel.pass-filtered.table.csv`

## 7. _**\[MODIFIED]**_**&#x20;Adjusted VAF with IMPACT FACETS**

* Uses FACETS copy number to compute adjusted VAF for **clonal variants**.
* Assumes normal copy number = 2 in all patients, and for chromosome X in females, But normal copy number = 1 for chromosomes X and Y in males.
* _**\[NEW]**_ Creates two outputs (only contain variants that overlap with the FACETS CCF MAF):
  1. Adjusted VAFs using all IMPACT samples. `{patient_id}-SNV-INDEL.allele_counts.hotspot_ch.filter.adj_vaf_all_impact.csv`
  2. Adjusted VAF from a single IMPACT sample selected for most variant overlap with ACCESS tumors.`{patient_id}.snv_indel.pass-filtered.adj_vaf.csv`

## Output files

Output files are organized as:

```
{outdir}/
└── intermediary/
│   └── small_variants/
│       └── {patient_id}/
│           ├── {patient_id}_all_small_calls.maf
│           ├── {patient_id}_genotyping_input.tsv
│           ├── {patient_id}-SNV-INDEL.allele_counts.csv
│           ├── {patient_id}-SNV-INDEL.allele_counts.hotspot_ch.csv
│           ├── {patient_id}-SNV-INDEL.allele_counts.hotspot_ch.filter.csv
│           ├── {patient_id}-SNV-INDEL.allele_counts.hotspot_ch.filter.adj_vaf_all_impact.csv
│           └── genotyped_mafs/
└── final/
    └── {patient_id}/
        ├── {patient_id}.snv_indel.pass-filtered.csv
        ├── {patient_id}.snv_indel.pass-filtered.table.csv
        └── {patient_id}.snv_indel.pass-filtered.adj_vaf.csv
```
