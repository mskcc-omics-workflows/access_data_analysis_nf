process MSI {
    tag "$patient_id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"

    input:
    tuple path(patient_json), val(patient_id)
    val reseach_access_msi_template
    path clinical_access_msi_file
    path clinical_impact_msi_file

    publishDir "${params.outdir}/final/${patient_id}", mode: 'copy', pattern: '*msi.csv'

    output:
        tuple path(patient_json), path("*msi.csv"), emit: msi_results

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ${workflow.projectDir}/bin/msi_analysis.py \\
        --patient_json $patient_json \\
        --reseach_access_msi_template $reseach_access_msi_template \\
        --clinical_access_msi_file $clinical_access_msi_file \\
        --clinical_impact_msi_file $clinical_impact_msi_file \\
        --output ${patient_id}.msi.csv
    """

}
