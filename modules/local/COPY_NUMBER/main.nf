process COPY_NUMBER {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"

    input:
    path patient_json
    val research_access_cna_template
    path clinical_cna_file

    publishDir "${params.outdir}/final_results/copy_number_variants", mode: 'copy', pattern: '*_CNA_calls.csv'

    output:
        tuple path(patient_json), path('*_CNA_calls.csv'), emit: copy_number_results

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ../../../bin/get_cna_calls.py \\
        --patient_json $patient_json \\
        --research_access_cna_template $research_access_cna_template \\
        --clinical_cna_file $clinical_cna_file \\
    """

}
