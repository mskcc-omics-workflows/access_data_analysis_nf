<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-accessanalysis_logo_dark.png">
    <img alt="nf-core/accessanalysis" src="docs/images/nf-core-accessanalysis_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/accessanalysis/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/accessanalysis/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/accessanalysis/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/accessanalysis/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/accessanalysis/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/accessanalysis)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23accessanalysis-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/accessanalysis)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

# ACCESS Data Analysis Pipeline

## Overview
Pipeline for analyzing variants from research and clinical ACCESS and IMPACT samples, and creating per-patient summaries of variants across the different samples.

## Dependencies
- Nextflow ≥ 24.04.2
- Python ≥ 3.9
- pandas
- numpy

## Usage

```bash
nextflow run main.nf -c nextflow.config -profile conda,juno,accessv1
```
- ***Note 1:** `juno` and `accessv1` are set as defaults and may be omitted.*
- ***Note 2: [MODIFIED]** The nextflow.config excludes JUNO-specific file path parameters, which are contained in `conf/juno.config`.*

## Input and Output parameters

*Set these in the `nextflow.config` file or specify on the command line:*
- **input** *(required)*: `/path/to/input/patient_id_mapping.csv`  
  CSV file with columns: `cmo_patient_id`, `dmp_patient_id`, `sex` of patients to include.
- **keep_research_samples_file** *(optional)*: File containing **CMO sample IDs** to include (no header).
  - If provided, other CMO samples for the patients will be ignored.
  - Use this to provide the same sample list you would specify in the R pipeline masterfile.
  - Use `null` to not restrict CMO samples to a list.
- **exclude_samples_file** *(optional)*: File containing **CMO and DMP sample IDs** to exclude (no header).
  - Exclusions override inclusions for CMO sample IDs.
  - Use `null` to not exclude any samples.
- **outdir** *(required)*: Output directory.

## Pipeline Steps

### 1. Research and Clinical ACCESS and IMPACT sample inference
- Parses ACCESS research BAM directories and key files for ACCESS/IMPACT clinical samples.  
- Uses include and exclude lists provided by user to keep or exclude samples. 
- Sets combined patient ID as `{cmo_patient_id}_{dmp_patient_id}`
  - [or `cmo_patient_id` or `dmp_patient_id` if other ID not present]
- Output files: `{patient_id}.json`.

**The subsequent steps are run separately for each patient `{patient_id}`.**
-----------------------------
### 2. Find FACETS fit
- Selects optimal FACETS fit for each IMPACT sample (prioritizing reviewed/QC-passed).  
- Falls back to default fit if needed.  
- Outputs tab-separated mapping of samples to FACETS fits.  
- Marks missing cases with `"MISSING"`.
- Output file: `{patient_id}_facets_fit.txt`
  
### 3. SNV/INDEL analysis steps

#### 3.1 Union MAF Generation
- Collects variants from ACCESS & IMPACT samples.  
- Produces unified MAF file with all variants.  
- ***[NEW]*** Does not filter any variants at this point 
- ***[NEW]*** Adds columns:  
  - **Called_In**: Samples containing the variant (semicolon-separated).  
  - **Clinical**: Marks variants from clinical samples as `"Signed Out"`
- Output file: `{patient_id}_all_small_calls.maf`

#### 3.2 Input generation for genotype variants
- Builds genotyping input with inferred BAM paths:  
  - Duplex & Simplex for ACCESS tumors.  
  - Unfiltered for ACCESS normals.  
  - Standard for IMPACT.
- Uses union MAF from previous step
- Output file: `{patient_id}_genotyping_input.tsv`
    
#### 3.3 Genotype Variants
- Runs `genotype_variants small_variants multiple-samples`.  
- Deduplication enabled.  
- Requires `GetBaseCountsMultiSample`.
- Output files: `{patient_id}/genotyped_mafs/*.maf`

#### 3.4 Aggregate variant allele counts across samples
- Produces per-sample allele counts + VAF for each variant.
- If there are M variants and N samples, there will be MxN rows in this table.
- Uses fragment counts for alt and total counts for ACCESS tumor samples and raw counts for ACCESS normal and IMPACT samples.
- Fragment counts are from the following BAM files:
  - ACCESS tumors: SIMPLEX + DUPLEX
  - ACCESS normals: UNFILTERED
  - IMPACT: STANDARD
- ***[NEW]*** Assigns call status per sample:
  - "called": variant was called in that sample
  - "low_coverage": total count is below coverage threshold (default - ACCESS: 100x, IMPACT: 50x)
  - "genotyped": previous conditions not met and VAF > 0
  - "": (empty) none of the above conditions met
- Output file: `{patient_id}-SNV-INDEL.allele_counts.csv`
      
#### 3.5 Annotate with Hotspot and CH lists
- Annotates against hotspot and clonal hematopoiesis lists.  
- Adds **Hotspot** and **CH** columns ("yes" if matched, empty if not).
- Output file: `{patient_id}-SNV-INDEL.allele_counts.hotspot_ch.csv`

#### 3.6 Filter variants ***[MODIFIED]***
- ***[NEW]*** Adds 'filter' column to variant table, marking variants for exclusion based on:
  - Presence in excluded genes list ('excluded_gene')
  - Matching excluded variant classifications ('excluded_classification')
  - ***[NEW]*** Low coverage (low_coverage) in ALL (tumor and normal, research and clinical) ACCESS samples ('low_access_cov')
  - ***[MODIFIED]*** For non-signed out variants if none of the ACCESS tumor samples meet minimum duplex alt threshold (default threshold: 3 for hotspot mutation, 5 for non-hotspot mutations) ('max_duplex_alt_count')
  - ***[MODIFIED]*** For non-signed out variants if the ratio of max VAF of tumor samples to max VAF of normal samples if less than threshold (default: 2) ('low_tumor_to_normal_vaf_ratio')
  - Multiple filter reasons are combined with semicolons (e.g., "excluded_gene;low_access_cov")
- Creates three output files:
  1. ***[NEW]*** Full variant list with annotations. `{patient_id}-SNV-INDEL.allele_counts.hotspot_ch.filter.csv`
  2. Filtered PASS-only list. `{patient_id}.snv_indel.pass-filtered.csv`
  3. A table format with one row per variant of the pass-filtered file. `{patient_id}.snv_indel.pass-filtered.table.csv`

#### 3.7 ***[MODIFIED]*** Adjusted VAF with IMPACT FACETS
- Uses FACETS copy number to compute adjusted VAF for **clonal variants**.
- ***[FIXED]*** Assumes normal copy number = 2 in all patients, and for chromosome X in females, But normal copy number = 1 for chromosomes X and Y in males.
- ***[NEW]*** Creates two outputs (only contain variants that overlap with the FACETS CCF MAF):
  1. Adjusted VAFs using all IMPACT samples. `{patient_id}-SNV-INDEL.allele_counts.hotspot_ch.filter.adj_vaf_all_impact.csv`
  2. Adjusted VAF from a single IMPACT sample selected for most variant overlap with ACCESS tumors.`{patient_id}.snv_indel.pass-filtered.adj_vaf.csv`

## Output Files

Results are organized as:
```
{outdir}/
├── intermediary/
│   └──facets_fit/
│           ├── {patient_id}_facets_fit.txt
│   └──patient_JSONs/
│           ├── {patient_id}_all_samples.json
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
        └── {patient_id}.snv_indel.pass-filtered.csv
        └── {patient_id}.snv_indel.pass-filtered.table.csv
        └── {patient_id}.snv_indel.pass-filtered.adj_vaf.csv
```


## Credits
- **@shguturu** - Original author
- **@kanika-arora** - Supervision and modifications

We thank the following people for their extensive assistance in the development of this pipeline: **@buehlere**


<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#accessanalysis` channel](https://nfcore.slack.com/channels/accessanalysis) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/accessanalysis for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
