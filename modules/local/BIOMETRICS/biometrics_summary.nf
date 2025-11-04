process BIOMETRICS_SUMMARY {
    tag "$patient_id"
    label 'biometrics_summary'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/biometrics:0.2.16':
        'ghcr.io/msk-access/biometrics:0.2.16' }"

    input:
    tuple val(patient_id), path(biometrics_genotype_csv), path(biometrics_sexmismatch_csv)

    output:
    tuple val(patient_id), path("${patient_id}.biometrics_summary.csv"), emit: biometrics_summary
    // Publish the summary CSV to the final results directory
    publishDir "${params.outdir}/final/${patient_id}", mode: 'copy', pattern: '*.biometrics_summary.csv'

    script:
    """
    python ${workflow.projectDir}/bin/summarize_biometrics_results.py \\
        --genotype ${biometrics_genotype_csv} \\
        --sexmismatch ${biometrics_sexmismatch_csv} \\
        --output ${patient_id}.biometrics_summary.csv
    """
}

