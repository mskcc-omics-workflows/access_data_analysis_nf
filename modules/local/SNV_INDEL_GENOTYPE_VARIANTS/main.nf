process SNV_INDEL_GENOTYPE_VARIANTS {
    tag "$patient_id"
    label 'genotype_variants'
    errorStrategy 'terminate'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/genotype_variants:sha-f0411c85':
        'ghcr.io/msk-access/genotype_variants:sha-f0411c85' }"

    input:
    tuple path(patient_json), val(patient_id), val(genotyping_input)
    val fasta_ref

    publishDir "${params.outdir}/intermediate/small_variants/${patient_id}/genotyped_mafs", mode: 'copy', pattern: '*.maf'

    output:
        tuple path(patient_json), val(patient_id), path("*.maf"), emit: genotyped_mafs
        stdout

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    
    genotype_variants small_variants multiple-samples \\
    -i ${genotyping_input} \\
    -r ${fasta_ref} \\
    --filter-duplicate 1 \\
    -g /juno/work/access/production/resources/tools/GetBaseCountsMultiSample/current/GetBaseCountsMultiSample \\
    -t ${task.cpus} \\

    """

}
