process SNV_INDEL_CREATE_GENOTYPE_INPUT {
    tag "$patient_id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6':
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6' }"

    input:
    tuple path(patient_sheet), val(patient_id), path(all_calls_maf)

    publishDir "${params.outdir}/intermediate/${patient_id}", mode: 'copy', pattern: '*genotyping_input.tsv'

    output:
        tuple path(patient_sheet), val(patient_id), path("*genotyping_input.tsv"), emit: genotyping_input

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${workflow.projectDir}/bin/genotype_variants_input.py \\
        --patient_sheet $patient_sheet \\
        --all_calls_maf $all_calls_maf
    """

}
