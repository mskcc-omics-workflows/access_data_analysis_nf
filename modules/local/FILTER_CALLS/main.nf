process FILTER_CALLS {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"

    input:
    path patient_json
    path mafs
    path facets_file

    publishDir 'output/final_results/small_variants', mode: 'copy'

    output:
        path "*.csv", emit: snv_table


    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ../../../bin/filter_calls.py \\
        --patient_json $patient_json \\
        --facets_fit $facets_file \\

    """

}
