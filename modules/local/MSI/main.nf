process MSI {
    tag "$patient_id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6':
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6' }"

    input:
    tuple path(patient_sheet), val(patient_id)
    path clinical_access_msi_file
    path clinical_impact_msi_file

    publishDir "${params.outdir}/final/${patient_id}", mode: 'copy', pattern: '*msi.csv'

    output:
        tuple path(patient_sheet), path("*msi.csv"), emit: msi_results

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ${workflow.projectDir}/bin/msi_analysis.py \\
        --patient_sheet $patient_sheet \\
        --clinical_access_msi_file $clinical_access_msi_file \\
        --clinical_impact_msi_file $clinical_impact_msi_file \\
        --output ${patient_id}.msi.csv
    """

}
