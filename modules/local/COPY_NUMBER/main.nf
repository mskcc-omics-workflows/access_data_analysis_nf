process COPY_NUMBER {
    tag "$patient_id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6':
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6' }"

    input:
    tuple path(patient_sheet), val(patient_id)
    path clinical_cna_file
    val access_copy_number_gene_list_v1
    val access_copy_number_gene_list_v2
    val p_value_threshold


    publishDir "${params.outdir}/intermediate/${patient_id}", mode: 'copy', pattern: "*cnv.csv"
    publishDir "${params.outdir}/final/${patient_id}", mode: 'copy', pattern: "*cnv.pass-filtered.csv"

    output:
        tuple path(patient_sheet), path('*cnv*.csv'), emit: copy_number_results

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ${workflow.projectDir}/bin/copy_number_variant_analysis.py \\
        --patient_sheet $patient_sheet \\
        --clinical_cna_file $clinical_cna_file \\
        --access_copy_number_gene_list_v1 $access_copy_number_gene_list_v1 \\
        --access_copy_number_gene_list_v2 $access_copy_number_gene_list_v2 \\
        --p_value_threshold $p_value_threshold \\
        --output ${patient_id}.cnv.csv \\
        --output_final ${patient_id}.cnv.pass-filtered.csv
    """

}
