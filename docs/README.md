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