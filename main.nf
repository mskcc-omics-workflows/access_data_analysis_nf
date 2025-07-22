#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/accessanalysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/accessanalysis
    Website: https://nf-co.re/accessanalysis
    Slack  : https://nfcore.slack.com/channels/accessanalysis
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ACCESSANALYSIS  } from './workflows/accessanalysis'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_accessanalysis_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_accessanalysis_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_accessanalysis_pipeline'
include { INFER_SAMPLES         } from './modules/local/INFER_SAMPLES/main'
include { GENOTYPE_VARIANTS_INPUT         } from './modules/local/GENOTYPE_VARIANTS_INPUT/main'
include { GENOTYPE_VARIANTS         } from './modules/local/GENOTYPE_VARIANTS/main'
include { GENERATE_MAF         } from './modules/local/GENERATE_MAF/main'
include { FIND_FACETS_FIT         } from './modules/local/FIND_FACETS_FIT/main'
include { FILTER_CALLS         } from './modules/local/FILTER_CALLS/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow MSK_ACCESS_DATA_ANALYSIS_NF {

    //take:
    take:
    patient_json
    //samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //

    GENERATE_MAF(
        patient_json,
        params.maf_template,
        params.dmp_calls_path
    )

    GENOTYPE_VARIANTS_INPUT(
        patient_json,
        params.fasta_ref,
        params.fasta_index,
        GENERATE_MAF.out.maf_results,

        // Research templates
        params.research_duplex_bam,
        params.research_duplex_bai,
        params.research_simplex_bam,
        params.research_simplex_bai,

        // Clinical access templates
        params.clinical_access_duplex_bam,
        params.clinical_duplex_bai,
        params.clinical_simplex_bam,
        params.clinical_simplex_bai,

        // Impact templates
        params.impact_standard_bam,
        params.impact_standard_bai
    )

    GENOTYPE_VARIANTS(
        patient_json,
        GENOTYPE_VARIANTS_INPUT.out.genotyping_input,
        params.fasta_ref,
        params.fasta_index
    )

    FIND_FACETS_FIT(
        params.facet_path,
        patient_json

    )

    FILTER_CALLS(
        patient_json,
        GENOTYPE_VARIANTS.out.genotyped_mafs,
        FIND_FACETS_FIT.out.facet_fit
    )

    //ACCESSANALYSIS (
    //    samplesheet
    //)
    emit:
    multiqc_report = null // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //

    INFER_SAMPLES (
        PIPELINE_INITIALISATION.out.samplesheet,
        params.include_samples_file,
        params.exclude_samples_file,
        params.dmp_access_key_path,
        params.dmp_impact_key_path,
        params.research_access_dir_template,
        params.access_sample_regex_pattern,
        params.impact_sample_regex_pattern
    )

    json_files = INFER_SAMPLES.out.all_samples_json.flatten()
    json_files | MSK_ACCESS_DATA_ANALYSIS_NF
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        MSK_ACCESS_DATA_ANALYSIS_NF.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
