process BIOMETRICS_QC_RESULT {
    tag "$patient_id"
    label 'biometrics_qc_result'

    input:
    tuple val(patient_id), path(biometrics_genotype_csv)

    output:
    tuple val(patient_id), path("${patient_id}.biometrics_qc_result.csv"), emit: biometrics_qc_result

    // Only publish the QC result CSV
    publishDir "${params.outdir}/final_results/biometrics", mode: 'copy', pattern: '*.biometrics_qc_result.csv'

    script:
    """
    qc_result=""

    # Determine QC status
    if grep -q "Unexpected" "${biometrics_genotype_csv}"; then
        qc_result="FAIL"
        echo "Patient \${patient_id} failed: Unexpected Mismatch detected in biometrics genotype comparison"
    elif [ "\$(wc -l < "${biometrics_genotype_csv}")" -le 1 ]; then
        qc_result="UNKNOWN"
        echo "Please check biometrics genotype comparison file for patient \${patient_id}"
    else
        qc_result="PASS"
    fi

    # Write QC result to file
    echo "patient_id,qc_result" > "${patient_id}.biometrics_qc_result.csv"
    echo "${patient_id},\${qc_result}" >> "${patient_id}.biometrics_qc_result.csv"
    """
}

