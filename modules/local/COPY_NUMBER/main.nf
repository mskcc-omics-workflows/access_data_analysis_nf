process COPY_NUMBER {
    tag "$patient_id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"

    input:
    tuple path(patient_json), val(patient_id)
    val research_access_cna_template
    path clinical_cna_file
    val access_copy_number_gene_list
    val p_value_threshold


    publishDir "${params.outdir}/intermediate/${patient_id}", mode: 'copy', pattern: "*cnv.csv"
    publishDir "${params.outdir}/final/${patient_id}", mode: 'copy', pattern: "*cnv.pass-filtered.csv"

    output:
        tuple path(patient_json), path('*cnv*.csv'), emit: copy_number_results

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ${workflow.projectDir}/bin/copy_number_variant_analysis.py \\
        --patient_json $patient_json \\
        --research_access_cna_template $research_access_cna_template \\
        --clinical_cna_file $clinical_cna_file \\
        --access_copy_number_gene_list $access_copy_number_gene_list \\
        --p_value_threshold $p_value_threshold \\
        --output ${patient_id}.cnv.csv \\
        --output_final ${patient_id}.cnv.pass-filtered.csv \\
    """

}
