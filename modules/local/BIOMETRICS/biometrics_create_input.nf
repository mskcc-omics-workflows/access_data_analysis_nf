process BIOMETRICS_CREATE_INPUT {
    tag "$patient_id"
    label 'process_single'

    input:
    tuple path(patient_json), val(patient_id)
    val research_access_standard_bam_template
    val clinical_access_standard_bam_template
    val clinical_impact_standard_bam_template

    publishDir "${params.outdir}/intermediate/biometrics/${patient_id}", mode: 'copy', pattern: '*biometrics_input.csv'

    output:
        tuple val(patient_id), path ("*biometrics_input.csv"), emit: biometrics_input

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${workflow.projectDir}/bin/create_biometrics_input_table.py \\
        --patient_json $patient_json \\
        --research_access_standard_bam_template $research_access_standard_bam_template \\
        --clinical_access_standard_bam_template $clinical_access_standard_bam_template \\
        --clinical_impact_standard_bam_template $clinical_impact_standard_bam_template \\

    """

}
