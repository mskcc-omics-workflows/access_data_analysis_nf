process BIOMETRICS_SEXMISMATCH {
    tag "$patient_id"
    label 'biometrics_genotype'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(patient_id), path(biometrics_input), path(biometrics_extract_db)

    output:
    tuple val(patient_id), path("${patient_id}.sex_mismatch.csv"), emit: biometrics_sexmismatch

    publishDir "${params.outdir}/intermediate/biometrics/${patient_id}", mode: 'copy', pattern: '*.sex_mismatch.csv'

    script:
    """

    # Run biometrics sex mismatch
    biometrics sexmismatch \\
        -i ${biometrics_input} \\
        -db ${biometrics_extract_db}

    mv sex_mismatch.csv ${patient_id}.sex_mismatch.csv
    """
}
