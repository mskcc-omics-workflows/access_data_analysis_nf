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
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_accessanalysis_pipeline'
include { getGenomeAttribute } from './subworkflows/local/utils_nfcore_accessanalysis_pipeline'
include { BIOMETRICS_CREATE_INPUT } from './modules/local/BIOMETRICS/biometrics_create_input'
include { BIOMETRICS_EXTRACT } from './modules/local/BIOMETRICS/biometrics_extract'
include { BIOMETRICS_GENOTYPE } from './modules/local/BIOMETRICS/biometrics_genotype'
include { BIOMETRICS_SEXMISMATCH } from './modules/local/BIOMETRICS/biometrics_sexmismatch'
include { BIOMETRICS_SUMMARY } from './modules/local/BIOMETRICS/biometrics_summary'
include { INFER_SAMPLES } from './modules/local/INFER_SAMPLES/main'
include { SNV_INDEL_CREATE_GENOTYPE_INPUT         } from './modules/local/SNV_INDEL_CREATE_GENOTYPE_INPUT/main'
include { SNV_INDEL_GENOTYPE_VARIANTS } from './modules/local/SNV_INDEL_GENOTYPE_VARIANTS/main'
include { SNV_INDEL_GENERATE_UNION_MAF } from './modules/local/SNV_INDEL_GENERATE_UNION_MAF/main'
include { FIND_FACETS_FIT } from './modules/local/FIND_FACETS_FIT/main'
include { SNV_INDEL_AGGREGATE_ALLELE_COUNTS } from './modules/local/SNV_INDEL_FILTER_CALLS/main'
include { SNV_INDEL_ADD_FACETS_ADJUSTED_VAF } from './modules/local/SNV_INDEL_FILTER_CALLS/main'
include { SNV_INDEL_ADD_FILTER_COL } from './modules/local/SNV_INDEL_FILTER_CALLS/main'
include { SNV_INDEL_ANNOTATE_HOTSPOT_CH } from './modules/local/SNV_INDEL_FILTER_CALLS/main'
include { COPY_NUMBER } from './modules/local/COPY_NUMBER/main'
include { STRUCTURAL_VARIANTS } from './modules/local/STRUCTURAL_VARIANTS/main'
include { MSI } from './modules/local/MSI/main'

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
 
    patient_meta = patient_json.map { file ->
        def json = new groovy.json.JsonSlurper().parseText(file.text)
        def patient_id = json.combined_id
        tuple(file, patient_id)
    }

    id_sex_mapping = patient_json.map { file ->
        def json = new groovy.json.JsonSlurper().parseText(file.text)
        def patient_id = json.combined_id
        def sex        = json.sex
    tuple(patient_id, sex)
}

    BIOMETRICS_CREATE_INPUT(
        patient_meta,
        params.file_paths.research_access.bam_file_template.standard,
        params.file_paths.clinical_access.bam_file_template.standard,
        params.file_paths.clinical_impact.bam_file_template.standard
    )
    
    BIOMETRICS_EXTRACT(
        BIOMETRICS_CREATE_INPUT.out.biometrics_input,
        params.fasta_ref,
        params.biometrics.bed,
        params.biometrics.vcf
    )
   
    BIOMETRICS_GENOTYPE(BIOMETRICS_EXTRACT.out.biometrics_extract)

    BIOMETRICS_SEXMISMATCH(BIOMETRICS_EXTRACT.out.biometrics_extract)

    BIOMETRICS_SUMMARY(
        BIOMETRICS_GENOTYPE.out.biometrics_genotype
        .join(BIOMETRICS_SEXMISMATCH.out.biometrics_sexmismatch, by:0)
    )

    SNV_INDEL_GENERATE_UNION_MAF(
        patient_meta,
        params.file_paths.research_access.variant_file_template.mutations,
        params.file_paths.clinical_impact.variant_file.mutations
    )


    SNV_INDEL_CREATE_GENOTYPE_INPUT(
        SNV_INDEL_GENERATE_UNION_MAF.out.maf_results,

        // Research ACCESS templates
        params.file_paths.research_access.bam_file_template.duplex,
        params.file_paths.research_access.bam_file_template.simplex,
        params.file_paths.research_access.bam_file_template.unfilter,

        // Clinical ACCESS templates
        params.file_paths.clinical_access.bam_file_template.duplex,
        params.file_paths.clinical_access.bam_file_template.simplex,
        params.file_paths.clinical_access.bam_file_template.unfilter,

        // Clinical IMPACT templates
        params.file_paths.clinical_impact.bam_file_template.standard
    )

    SNV_INDEL_GENOTYPE_VARIANTS(
        SNV_INDEL_CREATE_GENOTYPE_INPUT.out.genotyping_input,
        params.fasta_ref
    )

    FIND_FACETS_FIT(
        params.base_dirs.clinical_impact.facets_dir,
        patient_meta
    )
                            

    snv_indel_aggregate_input = SNV_INDEL_GENOTYPE_VARIANTS.out.genotyped_mafs
    .join(SNV_INDEL_GENERATE_UNION_MAF.out.maf_results, by: 1)
    .map { patient_id, geno_json, geno_mafs, union_json, union_maf -> 
        tuple(geno_json, patient_id, geno_mafs, union_maf)
    }
    
    SNV_INDEL_AGGREGATE_ALLELE_COUNTS(
        snv_indel_aggregate_input,
        params.variant_filter_rules.access_min_cov,
        params.variant_filter_rules.impact_min_cov
    )
    SNV_INDEL_ANNOTATE_HOTSPOT_CH(
        SNV_INDEL_AGGREGATE_ALLELE_COUNTS.out.snv_indel_results,
        params.hotspot_list,
        params.ch_list
    )
    SNV_INDEL_ADD_FILTER_COL(
        SNV_INDEL_ANNOTATE_HOTSPOT_CH.out.hotspot_ch_annotated_snv_indel, 
        params.variant_filter_rules.exclude_genes,
        params.variant_filter_rules.exclude_classifications,
        params.variant_filter_rules.hotspot_cutoff,
        params.variant_filter_rules.non_hotspot_cutoff,
        params.variant_filter_rules.vaf_ratio_threshold
    )
    
//    snv_indel_adj_vaf_input = SNV_INDEL_ADD_FILTER_COL.out.filtered_snv_indel
//    .join(FIND_FACETS_FIT.out.facets_fit, by: 0)
//    .map { patient_id, filtered_snv_indel, facets_fit -> 
//        tuple(patient_id, facets_fit, filtered_snv_indel)
//    }


    snv_indel_adj_vaf_input = id_sex_mapping
        .join(FIND_FACETS_FIT.out.facets_fit, by: 0)
        .join(SNV_INDEL_ADD_FILTER_COL.out.filtered_snv_indel, by: 0)
        .map { patient_id, sex, facets_fit, filtered_snv_indel ->
            tuple(patient_id, sex, facets_fit, filtered_snv_indel)
    }

    SNV_INDEL_ADD_FACETS_ADJUSTED_VAF(
        snv_indel_adj_vaf_input
    )


//    STRUCTURAL_VARIANTS(
//        patient_json,
//        params.file_paths.research_access.variant_file_template.sv,
//        params.file_paths.clinical_access.variant_file.sv,
//        params.file_paths.clinical_impact.variant_file.sv,
//        params.access_structural_variant_gene_list
//    )

//    MSI(
//        patient_json,
//        params.file_paths.research_access.variant_file_template.msi,
//        params.file_paths.clinical_access.variant_file.msi,
//        params.file_paths.clinical_impact.variant_file.msi
//    )

//    COPY_NUMBER (
//        patient_json,
//        params.file_paths.research_access.variant_file_template.cna,
//        params.file_paths.clinical_impact.variant_file.cna,
//        params.access_copy_number_gene_list,
//        params.research_access_copy_number_p_value_filter
//    )

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
        params.keep_research_samples_file ?: "$projectDir/assets/NO_INCLUDE_FILE",
        params.exclude_samples_file ?: "$projectDir/assets/NO_EXCLUDE_FILE",
        params.file_paths.clinical_access.key_file,
        params.file_paths.clinical_impact.key_file,
        params.base_dirs.research_access.bam_dir_template,
        params.clinical_access_sample_regex_pattern,
        params.clinical_impact_sample_regex_pattern
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
