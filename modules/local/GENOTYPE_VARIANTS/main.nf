process GENOTYPE_VARIANTS {
    tag "$patient_json"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/genotype_variants:0.3.9':
        'ghcr.io/msk-access/genotype_variants:0.3.9' }"

    input:
    path patient_json
    path genotyping_input
    path fasta_ref
    path fasta_index

    publishDir 'output/intermediary/genotyped_mafs', mode: 'copy'

    output:
    path "*.maf", emit: genotyped_mafs
    stdout

    when:
    task.ext.when == null || task.ext.when

    script:

    """

    genotype_variants small_variants multiple-samples \\
    -i ${genotyping_input} \\
    -r ${fasta_ref} \\
    --filter-duplicate 1 \\
    -g /work/access/production/resources/tools/GetBaseCountsMultiSample/current/GetBaseCountsMultiSample \\
    -t ${task.cpus} \\

    """

}
