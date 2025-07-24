process GENOTYPE_VARIANTS_INPUT {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"

    input:
    path patient_json
    path all_calls_maf
    val research_access_duplex_bam
    val research_access_simplex_bam
    val clinical_access_duplex_bam
    val clinical_access_simplex_bam
    val clinical_access_standard_bam
    val clinical_impact_standard_bam

    publishDir "${params.outdir}/intermediary/genotyping_input", mode: 'copy'

    output:
    path "*genotyping_input.tsv", emit: genotyping_input

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ../../../bin/genotype_variants_input.py \\
        --patient_json $patient_json \\
        --all_calls_maf $all_calls_maf \\
        --research_access_duplex_bam $research_access_duplex_bam \\
        --research_access_simplex_bam $research_access_simplex_bam \\
        --clinical_access_duplex_bam $clinical_access_duplex_bam \\
        --clinical_access_simplex_bam $clinical_access_simplex_bam \\
        --clinical_access_standard_bam $clinical_access_standard_bam \\
        --clinical_impact_standard_bam $clinical_impact_standard_bam \\

    """

}
