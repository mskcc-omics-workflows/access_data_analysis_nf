process FIND_FACETS_FIT {
    tag "$patient_json"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/genotype_variants:0.3.9':
        'ghcr.io/msk-access/genotype_variants:0.3.9' }"

    input:
    val facets_dir
    path patient_json

    publishDir "${params.outdir}/intermediary/facets_fit", mode: 'copy', pattern: '*facets_fit.txt'

    output:
        tuple path(patient_json), path("*facets_fit.txt"), emit: facets_fit

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ../../../bin/facets_fit.py \\
        --facets_dir $facets_dir \\
        --patient_json $patient_json \\
    """


}
