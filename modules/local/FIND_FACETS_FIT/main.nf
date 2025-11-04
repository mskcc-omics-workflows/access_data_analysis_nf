process FIND_FACETS_FIT {
    tag "$patient_id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6':
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6' }"

    input:
    val facets_dir
    tuple path(patient_json), val(patient_id)

    publishDir "${params.outdir}/intermediate/${patient_id}", mode: 'copy', pattern: '*facets_fit.txt'

    output:
        tuple val(patient_id), path("*facets_fit.txt"), emit: facets_fit

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ${workflow.projectDir}/bin/facets_fit.py \\
        --facets_dir $facets_dir \\
        --patient_json $patient_json \\
    """


}
