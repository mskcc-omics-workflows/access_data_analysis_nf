process BIOMETRICS_QC_RESULT {
    tag "$patient_id"
    label 'biometrics_qc_result'

    input:
    tuple path(patient_json), val(patient_id), path(genotype_csv)

    output:
    tuple path(patient_json), val(patient_id), val(biometrics_qc_result), path("${patient_id}.biometrics_qc_result.csv"), emit: biometrics_qc_result

    publishDir "${params.outdir}/final_results/biometrics", mode: 'copy', pattern: '*.biometrics_qc_result.csv' 
    script:
    """

    if grep -q "Unexpected" ${genotype_csv}; then
        biometrics_qc_result="FAIL"
        echo "Patient ${patient_id} failed: Unexpected Mismatch detected in biometrics genotype comparison"
    else
        biometrics_qc_result="PASS"
    fi
    echo "patient_id,qc_result" > ${patient_id}.biometrics_qc_result.csv
    echo "${patient_id},${biometrics_qc_result}" >> ${patient_id}.biometrics_qc_result.csv
    """
}
