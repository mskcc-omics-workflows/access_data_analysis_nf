process GENERATE_MAF {
    label 'process_single'

    input:
    path patient_json
    val research_access_mutations_maf_template
    path dmp_mutations_file
    val exclude_genes
    val exclude_classifications

    publishDir "${params.outdir}/intermediate/MAFs", mode: 'copy', pattern: '*_all_small_calls.maf'

    output:
        tuple path(patient_json), path("*_all_small_calls.maf"), emit: maf_results

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ../../../bin/generate_maf.py \\
        --patient_json $patient_json \\
        --research_access_mutations_maf_template $research_access_mutations_maf_template \\
        --dmp_mutations_file $dmp_mutations_file \\
        --exclude_genes $exclude_genes \\
        --exclude_classifications $exclude_classifications \\
    """

}
