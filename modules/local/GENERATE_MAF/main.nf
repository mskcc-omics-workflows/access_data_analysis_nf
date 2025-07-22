process GENERATE_MAF {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"

    input:
    path patient_json
    val maf_template
    path dmp_calls_path

    publishDir 'output/intermediary/MAFs', mode: 'copy'

    output:
        path "*_all_calls.maf", emit: maf_results


    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ../../../bin/generate_maf.py \\
        --patient_json $patient_json \\
        --maf_template $maf_template \\
        --dmp_calls_path $dmp_calls_path \\
    """

}
