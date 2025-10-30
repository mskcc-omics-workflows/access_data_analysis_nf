process SNV_INDEL_GENERATE_UNION_MAF {
    tag "$patient_id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6':
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.6' }"

    input:
    tuple path(patient_json), val(patient_id)
    val research_access_mutations_maf_template
    path dmp_mutations_file

    publishDir "${params.outdir}/intermediate/${patient_id}", mode: 'copy', pattern: '*snv_indel.union.maf'

    output:
        tuple path(patient_json), val(patient_id), path("*snv_indel.union.maf"), emit: maf_results

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ${workflow.projectDir}/bin/generate_snv_indel_union_maf.py \\
        --patient_json $patient_json \\
        --research_access_mutations_maf_template $research_access_mutations_maf_template \\
        --dmp_mutations_file $dmp_mutations_file \\
        --output ${patient_id}.snv_indel.union.maf \\
    """

}
