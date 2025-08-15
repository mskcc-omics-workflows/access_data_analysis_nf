process STRUCTURAL_VARIANTS {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"

    input:
    path patient_json
    val research_access_sv_template
    path clinical_access_sv_file, name: "clinical_access_data.txt"
    path clinical_impact_sv_file, name: "clinical_impact_data.txt"
    val access_structural_variant_gene_list

    publishDir "${params.outdir}/final_results/structural_variants", mode: 'copy', pattern: "*SV.csv"

    output:
        tuple path(patient_json), path("*SV.csv"), emit: sv_results

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ../../../bin/structural_variant_calls.py \\
        --patient_json $patient_json \\
        --research_access_sv_template $research_access_sv_template \\
        --clinical_access_sv_file $clinical_access_sv_file \\
        --clinical_impact_sv_file $clinical_impact_sv_file \\
        --access_structural_variant_gene_list $access_structural_variant_gene_list \\
    """

}
