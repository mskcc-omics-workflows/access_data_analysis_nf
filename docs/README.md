# nf-core/accessanalysis: Documentation

The nf-core/accessanalysis documentation is split into the following pages:

- [Usage](usage.md)
  - An overview of how the pipeline works, how to run it and a description of all of the different command-line flags.
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.

You can find a lot more documentation about installing, configuring and running nf-core pipelines on the website: [https://nf-co.re](https://nf-co.re)


## Input

The pipeline is driven by a **Voyager-prepared cohort samplesheet** (`--input`, a CSV). One
row per sample, with every BAM / MAF / variant-file path already resolved by Voyager. See
[`assets/schema_input.json`](../assets/schema_input.json) for the full column list and
[`assets/samplesheet.csv`](../assets/samplesheet.csv) for an example. `SPLIT_SAMPLESHEET`
slices this into one CSV per `combined_id` and the pipeline fans out one run per patient.

## Configuration

Config files in `/conf` supply only **cohort-level** resources — the reference FASTA, biometrics
BED/VCF, hotspot / CH lists, the FACETS directory, and the DMP aggregate files
(`data_mutations_extended.txt`, `data_CNA.txt`, `data_sv.txt`, `data_clinical_sample.txt`,
ADMIE MSI). Per-sample paths are **not** here — they come from the samplesheet.

### Profiles

- `voyager.config`: default; MSK Voyager environment
- `iris.config`: MSK Iris cluster
- `juno.config`: Juno/Terra cluster

### Running the Pipeline

```bash
nextflow run main.nf -profile voyager,singularity --input <cohort_samplesheet.csv> --outdir <outdir>
```

## Union MAF Generation
The script `generate_snv_indel_union_maf.py` aggregates and filters mutation calls from the
per-sample research MAFs (the `maf` column of the samplesheet) and the clinical DMP
`data_mutations_extended.txt`.

### Mutation File Format
Each mutation file (research or clinical) is expected to be tab-delimited and contain the following fields:
- `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`
- `Reference_Allele`, `Tumor_Seq_Allele1`, `Tumor_Seq_Allele2`
- `Tumor_Sample_Barcode`, `Mutation_Status`, `Status`
- `Variant_Classification`, `HGVSp`, `HGVSp_Short`

If file structure changes:
- `parse_mutation_file()` may fail to parse or filter mutations.

### Filtering Logic

During Parsing:
- Excludes **germline mutations**: `Mutation_Status == 'GERMLINE'`
- For **research mutations**: excludes any mutation where `Status != ''`
- For **clinical mutations**: includes only rows where `Tumor_Sample_Barcode` contains the `dmp_id` and skips metadata lines with `"sequenced_samples:"`

After Merging:
1. **Deduplication**:
   Removes duplicates based on:
   - `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Variant_Classification`, `Reference_Allele`, `Tumor_Seq_Allele2`

2. **Gene Filtering**:
   Removes any row where `Hugo_Symbol` **contains any substring** from `exclude_genes`

3. **Classification Filtering**:
   Removes any row where `Variant_Classification` **exactly matches** a string in `exclude_classifications`

Any changes in filtering should be made in the `nextflow.config`. Genes or classification types can be added to the exclude list as below, separated by commas.

``` bash
variant_filter_rules = [
    exclude_genes           : "RP11-,RET,BRAF"
    exclude_classifications : "Silent"
]
```

## Generating Input Table for Genotype Variants
The script `genotype_variants_input.py` creates the metadata TSV used by
`genotype_variants small_variants multiple-samples`. Each row is a sample to be genotyped
against the patient's union call set.

BAM paths are read directly from the samplesheet columns via `patient_sheet.py`:

- research / clinical ACCESS **tumor** → `duplex_bam` + `simplex_bam`
- research / clinical ACCESS **normal** → `standard_bam` := `unfilter_bam`
- clinical IMPACT → `standard_bam`

A sample is skipped (with a warning) if a BAM it needs is blank in the samplesheet.

### Output

A TSV file named `<combined_id>_genotyping_input.tsv` with columns:
- `patient_id` (combined patient id)
- `sample_id`
- `standard_bam`, `duplex_bam`, `simplex_bam`
- `maf` (full path to the patient union MAF)
