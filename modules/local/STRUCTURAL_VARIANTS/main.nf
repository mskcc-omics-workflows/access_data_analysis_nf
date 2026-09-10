process STRUCTURAL_VARIANTS {
    tag "$patient_id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6':
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6' }"

    input:
    tuple path(patient_sheet), val(patient_id)
    path clinical_impact_sv_file
    val access_structural_variant_gene_list

    publishDir "${params.outdir}/final/${patient_id}", mode: 'copy', pattern: "*sv.csv"

    output:
        tuple path(patient_sheet), path("*sv.csv"), emit: sv_results

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ${workflow.projectDir}/bin/structural_variant_analysis.py \\
        --patient_sheet $patient_sheet \\
        --clinical_sv_file $clinical_impact_sv_file \\
        --access_structural_variant_gene_list $access_structural_variant_gene_list \\
        --output ${patient_id}.sv.csv
    """

}
