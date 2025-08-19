process BIOMETRICS_EXTRACT {
    tag "$patient_id"
    label 'biometrics_extract'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(patient_id), path(biometrics_input)
    path fasta_ref
    path biometrics_bed
    path biometrics_vcf

    output:
    tuple val(patient_id), path(biometrics_input), path("${patient_id}"), emit: biometrics_extract

    publishDir "${params.outdir}/intermediary/biometrics/extract_db", mode: 'copy', pattern: "${patient_id}"

    script:
    """
    # write outputs there
    biometrics extract \
        -i ${biometrics_input} \
        -f ${fasta_ref} \
        --vcf ${biometrics_vcf} \
        --bed ${biometrics_bed} \
        -db ${patient_id} \
        -t ${task.cpus}
    """
}

