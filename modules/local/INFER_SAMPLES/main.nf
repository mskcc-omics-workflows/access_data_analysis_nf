process INFER_SAMPLES {
    label 'process_single'
    errorStrategy 'terminate'

    input:
    path id_mapping_file
    path keep_research_samples_file
    path exclude_samples_file
    path clinical_access_key_file, name: "access_key.txt"
    path clinical_impact_key_file, name: "dmp_key.txt"
    val research_access_bam_dir_template
    val clinical_access_sample_regex_pattern
    val clinical_impact_sample_regex_pattern        
    val research_access_manifest_file_template
    val research_access_mutations_maf_template

    publishDir "${params.outdir}/intermediate/patient_JSONs", mode: 'copy'

    output:
    path "*.json", emit: all_samples_json

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    command="python3 ../../../bin/infer_samples.py \\
        --id_mapping_file ${id_mapping_file} \\
        --clinical_access_key_file ${clinical_access_key_file} \\
        --clinical_impact_key_file ${clinical_impact_key_file} \\
        --research_access_bam_dir_template '${research_access_bam_dir_template}' \\
        --clinical_access_sample_regex_pattern '${clinical_access_sample_regex_pattern}' \\
        --clinical_impact_sample_regex_pattern '${clinical_impact_sample_regex_pattern}' \\
        --research_access_manifest_file_template '${research_access_manifest_file_template}' \\
        --research_access_mutations_maf_template '${research_access_mutations_maf_template}'"
    
    if [ -s "${keep_research_samples_file}" ]; then
        command="\$command --keep_research_samples_file '${keep_research_samples_file}'"
    fi
    
    if [ -s "${exclude_samples_file}" ]; then
        command="\$command --exclude_samples_file '${exclude_samples_file}'"
    fi
    
    echo \$command
    eval \$command
    """
}
