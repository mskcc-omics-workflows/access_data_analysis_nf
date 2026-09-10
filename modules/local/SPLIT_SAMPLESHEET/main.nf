process SPLIT_SAMPLESHEET {
    label 'process_single'
    errorStrategy 'terminate'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6':
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6' }"

    input:
    path samplesheet

    publishDir "${params.outdir}/intermediate/patient_samplesheets", mode: 'copy'

    output:
    path "*.samplesheet.csv", emit: patient_samplesheets

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${workflow.projectDir}/bin/split_samplesheet.py \\
        --samplesheet ${samplesheet}
    """
}
