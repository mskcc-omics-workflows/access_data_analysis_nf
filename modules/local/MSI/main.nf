process MSI {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"

    input:
    path patient_json
    val reseach_access_msi_template
    path clinical_access_msi_file
    path clinical_impact_msi_file

    publishDir "${params.outdir}/final_results/microsatellite_instability", mode: 'copy', pattern: '*_MSI.csv'

    output:
        tuple path(patient_json), path("*_MSI.csv"), emit: msi_results

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ../../../bin/msi_scores.py \\
        --patient_json $patient_json \\
        --reseach_access_msi_template $reseach_access_msi_template \\
        --clinical_access_msi_file $clinical_access_msi_file \\
        --clinical_impact_msi_file $clinical_impact_msi_file \\
    """

}
