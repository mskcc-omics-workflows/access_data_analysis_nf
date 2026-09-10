process BIOMETRICS_CREATE_INPUT {
    tag "$patient_id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6':
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6' }"

    input:
    tuple path(patient_sheet), val(patient_id)

    publishDir "${params.outdir}/intermediate/biometrics/${patient_id}", mode: 'copy', pattern: '*biometrics_input.csv'

    output:
        tuple val(patient_id), path ("*biometrics_input.csv"), emit: biometrics_input

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${workflow.projectDir}/bin/create_biometrics_input_table.py \\
        --patient_sheet $patient_sheet
    """

}
