process SNV_INDEL_CREATE_GENOTYPE_INPUT {
    tag "$patient_id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6':
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6' }"

    input:
    tuple path(patient_json), val(patient_id), path(all_calls_maf)
    val research_access_duplex_bam_template
    val research_access_simplex_bam_template
    val research_access_unfilter_bam_template
    val clinical_access_duplex_bam_template
    val clinical_access_simplex_bam_template
    val clinical_access_unfilter_bam_template
    val clinical_impact_standard_bam_template

    publishDir "${params.outdir}/intermediate/small_variants/${patient_id}", mode: 'copy', pattern: '*genotyping_input.tsv'

    output:
        tuple path(patient_json), val(patient_id), path("*genotyping_input.tsv"), emit: genotyping_input

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${workflow.projectDir}/bin/genotype_variants_input.py \\
        --patient_json $patient_json \\
        --all_calls_maf $all_calls_maf \\
        --research_access_duplex_bam_template $research_access_duplex_bam_template \\
        --research_access_simplex_bam_template $research_access_simplex_bam_template \\
        --research_access_unfilter_bam_template $research_access_unfilter_bam_template \\
        --clinical_access_duplex_bam_template $clinical_access_duplex_bam_template \\
        --clinical_access_simplex_bam_template $clinical_access_simplex_bam_template \\
        --clinical_access_unfilter_bam_template $clinical_access_unfilter_bam_template \\
        --clinical_impact_standard_bam_template $clinical_impact_standard_bam_template \\

    """

}
