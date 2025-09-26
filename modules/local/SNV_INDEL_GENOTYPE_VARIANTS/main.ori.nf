process GENOTYPE_VARIANTS {
    tag "$patient_json"
    label 'genotype_variants'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/genotype_variants:0.3.9':
        'ghcr.io/msk-access/genotype_variants:0.3.9' }"

    input:
    tuple path(patient_json), val(genotyping_input)
    val fasta_ref

    publishDir "${params.outdir}/intermediary/small_variants/genotyped_mafs", mode: 'copy', pattern: '*.maf'

    output:
        tuple path(patient_json), path("*.maf"), emit: genotyped_mafs
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
