process FIND_FACETS_FIT {
    tag "$patient_json"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/genotype_variants:0.3.9':
        'ghcr.io/msk-access/genotype_variants:0.3.9' }"

    input:
    val facet_path
    path patient_json

    publishDir 'output/intermediary/facets_fit', mode: 'copy'

    output:
    path "*facets_maf_path.txt", emit: maf_path

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ../../../bin/facets_fit.py \\
        --facet_path $facet_path \\
        --patient_json $patient_json \\
        --best_fit "True" \\
    """


}
