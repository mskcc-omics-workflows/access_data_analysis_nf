process SNV_INDEL_GENOTYPE_VARIANTS {
    tag "$patient_id"
    label 'genotype_variants'
    errorStrategy 'terminate'

    conda "${moduleDir}/environment.yml"
    // 0.3.12 fixes the builtin-all() shadowing in generate_gbcms_cmd that broke
    // `multiple-samples` on 0.3.6-0.3.11. Bundles GetBaseCountsMultiSample 1.2.5
    // at /usr/local/bin, so no external gbcms_path is needed.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/genotype_variants:0.3.12':
        'ghcr.io/msk-access/genotype_variants:0.3.12' }"

    input:
    tuple path(patient_sheet), val(patient_id), val(genotyping_input)
    val fasta_ref

    publishDir "${params.outdir}/intermediate/${patient_id}/genotyped_mafs", mode: 'copy', pattern: '*_genotyped.maf'

    output:
        tuple path(patient_sheet), val(patient_id), path("*.maf"), emit: genotyped_mafs
        stdout

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    genotype_variants small_variants multiple-samples \\
        -i ${genotyping_input} \\
        -r ${fasta_ref} \\
        -g /usr/local/bin/GetBaseCountsMultiSample \\
        --filter-duplicate 1 \\
        -t ${task.cpus}
    """

}
