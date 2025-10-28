process STRUCTURAL_VARIANTS {
    tag "$patient_id"
    label 'process_single'

    input:
    tuple path(patient_json), val(patient_id)
    val research_access_sv_template
    path clinical_impact_sv_file
    val access_structural_variant_gene_list

    publishDir "${params.outdir}/final/${patient_id}", mode: 'copy', pattern: "*sv.csv"

    output:
        tuple path(patient_json), path("*sv.csv"), emit: sv_results

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ${workflow.projectDir}/bin/structural_variant_analysis.py \\
        --patient_json $patient_json \\
        --research_access_sv_template $research_access_sv_template \\
        --clinical_sv_file $clinical_impact_sv_file \\
        --access_structural_variant_gene_list $access_structural_variant_gene_list \\
        --output ${patient_id}.sv.csv \\
    """

}
