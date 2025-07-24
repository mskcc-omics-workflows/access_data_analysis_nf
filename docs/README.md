# nf-core/accessanalysis: Documentation

The nf-core/accessanalysis documentation is split into the following pages:

- [Usage](usage.md)
  - An overview of how the pipeline works, how to run it and a description of all of the different command-line flags.
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.

You can find a lot more documentation about installing, configuring and running nf-core pipelines on the website: [https://nf-co.re](https://nf-co.re)


## Configuration

This pipeline uses configuration files to manage all path-specific settings. All paths to reference data, bam files, variant files, and resources are defined in the configuration files located in the `/conf` directory.

### Default Configs

- `juno.config`: Configuration for running in Juno/Terra cluster  
- `iris.config`: Configuration for running in Iris

Each config file defines:

- Input file paths (e.g., sample mappings, inclusion/exclusion lists)  
- Regex patterns for sample filtering  
- Directories and templates for locating BAM and variant files  
- Resource paths for reference genomes and annotation files  

You **should not** need to modify the pipeline script directly.

### Running the Pipeline

To run the pipeline with the Juno config, use:

```bash
nextflow run main.nf -c nextflow.config -process.echo -profile juno
```

## MAF Generation
The script `generate_maf.py` aggregates and filters mutation calls from both research and clinical maf files. 

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
The script `genotype_variants_input.py` creates the metadata table used by `genotype_variants`. 
Each row is a sample to be genotyped, and the columns are simplex, duplex, or standard bam paths associated with the sample.
The script `infer_bams.py` is a helper script that generates the bam paths given a template file and sample metadata. 

### Assumptions

- The research bam path templates use placeholders: `sample_id` and `cmo_patient_id` 
- The clincal bam path tempaltes use placeholders: `anon_id`, `anon_id_fl`, and `anon_id_sl` (referring to the first and second letter of the anon id)
- The base paths of the bams can be changed in the `nextflow.config` file, but the script will expect these placeholders to exist in the paths.

- BAM path templates:
  - `--research_access_duplex_bam`
  - `--research_access_simplex_bam`
  - `--clinical_access_duplex_bam`
  - `--clinical_access_simplex_bam`
  - `--clinical_access_standard_bam`
  - `--clinical_impact_standard_bam`


### Output

A TSV file named `<combined_id>_genotyping_input.tsv` with columns:
- `sample_id`
- `duplex_bam`, `simplex_bam`, or `standard_bam`
- `patient_id` (combined patient id)
- `maf` (full path to input MAF)

### BAM Path Validation

- The BAM path and index file path are validated in `infer_bams.py`.
- If the BAM file or its `.bai` index are missing, the path is replaced with `"MISSING_PATH"` in the input table and a warning is printed.
