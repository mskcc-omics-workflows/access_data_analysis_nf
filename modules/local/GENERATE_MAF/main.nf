process GENERATE_MAF {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"

    input:
    path patient_json
    val research_access_mutations_maf
    path dmp_mutations_file
    val exclude_genes
    val exclude_classifications

    publishDir "${params.outdir}/intermediary/MAFs", mode: 'copy'

    output:
        path "*_all_small_calls.maf", emit: maf_results


    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ../../../bin/generate_maf.py \\
        --patient_json $patient_json \\
        --research_access_mutations_maf $research_access_mutations_maf \\
        --dmp_mutations_file $dmp_mutations_file \\
        --exclude_genes $exclude_genes \\
        --exclude_classifications $exclude_classifications \\
    """

}
