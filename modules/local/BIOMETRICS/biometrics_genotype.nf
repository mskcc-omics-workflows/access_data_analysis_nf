process BIOMETRICS_GENOTYPE {
    tag "$patient_id"
    label 'biometrics_genotype'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/biometrics:0.2.16':
        'ghcr.io/msk-access/biometrics:0.2.16' }"

    input:
    tuple val(patient_id), path(biometrics_input), path(biometrics_extract_db)

    output:
    tuple val(patient_id), path("${patient_id}.genotype_comparison.csv"), emit: biometrics_genotype

    publishDir "${params.outdir}/intermediate/biometrics/${patient_id}", mode: 'copy', pattern: '*.genotype_comparison.csv'

    script:
    """

    # Run biometrics genotype in that directory
    biometrics genotype \\
        -i ${biometrics_input} \\
        -db ${biometrics_extract_db}

    mv genotype_comparison.csv ${patient_id}.genotype_comparison.csv
    """
}
