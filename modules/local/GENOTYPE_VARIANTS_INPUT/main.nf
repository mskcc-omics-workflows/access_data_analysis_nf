process GENOTYPE_VARIANTS_INPUT {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"

    input:
    path patient_json
    path fasta_ref
    path fasta_index
    path maf_output

    // all template paths
    val research_duplex_bam
    val research_duplex_bai
    val research_simplex_bam
    val research_simplex_bai

    val clinical_access_duplex_bam
    val clinical_duplex_bai
    val clinical_simplex_bam
    val clinical_simplex_bai

    val impact_standard_bam
    val impact_standard_bai

    publishDir 'output/intermediary/genotyping_input', mode: 'copy'

    output:
    path "*genotyping_input.tsv", emit: genotyping_input

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ../../../bin/genotype_variants_input.py \\
        --patient_json $patient_json \\
        --fasta_ref $fasta_ref \\
        --maf_results $maf_output \\
        --threads ${task.cpus} \\
        --research_duplex_bam $research_duplex_bam \\
        --research_duplex_bai $research_duplex_bai \\
        --research_simplex_bam $research_simplex_bam \\
        --research_simplex_bai $research_simplex_bai \\
        --clinical_access_duplex_bam $clinical_access_duplex_bam \\
        --clinical_duplex_bai $clinical_duplex_bai \\
        --clinical_simplex_bam $clinical_simplex_bam \\
        --clinical_simplex_bai $clinical_simplex_bai \\
        --impact_standard_bam $impact_standard_bam \\
        --impact_standard_bai $impact_standard_bai
    """

}
