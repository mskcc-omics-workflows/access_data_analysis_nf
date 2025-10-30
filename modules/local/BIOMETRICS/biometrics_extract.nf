process BIOMETRICS_EXTRACT {
    tag "$patient_id"
    label 'biometrics_extract'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/biometrics:0.2.16':
        'ghcr.io/msk-access/biometrics:0.2.16' }"

    input:
    tuple val(patient_id), path(biometrics_input)
    path fasta_ref
    path biometrics_bed
    path biometrics_vcf

    output:
    tuple val(patient_id), path(biometrics_input), path("extract_db"), emit: biometrics_extract

    publishDir "${params.outdir}/intermediate/biometrics/${patient_id}", mode: 'copy', pattern: "extract_db"

    script:
    """
    # write outputs there
    biometrics extract \
        -i ${biometrics_input} \
        -f ${fasta_ref} \
        --vcf ${biometrics_vcf} \
        --bed ${biometrics_bed} \
        -db 'extract_db' \
        -t ${task.cpus}
    """
}

