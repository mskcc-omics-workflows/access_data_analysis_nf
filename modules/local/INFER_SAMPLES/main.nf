process INFER_SAMPLES {
    label 'process_single'

    input:
    path id_mapping_file
    path keep_research_samples_file
    path exclude_samples_file
    path clinical_access_key_file, name: "access_key.txt"
    path clinical_impact_key_file, name: "dmp_key.txt"
    val research_access_bam_dir_template
    val clinical_access_sample_regex_pattern
    val clinical_impact_sample_regex_pattern        

    publishDir "${params.outdir}/intermediate/patient_JSONs", mode: 'copy'

    output:
    path "*.json", emit: all_samples_json

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${workflow.projectDir}/bin/infer_samples.py \
        --id_mapping_file ${id_mapping_file} \
        --clinical_access_key_file ${clinical_access_key_file} \
        --clinical_impact_key_file ${clinical_impact_key_file} \
        --research_access_bam_dir_template '${research_access_bam_dir_template}' \
        --clinical_access_sample_regex_pattern '${clinical_access_sample_regex_pattern}' \
        --clinical_impact_sample_regex_pattern '${clinical_impact_sample_regex_pattern}' \
        ${keep_research_samples_file ? "--keep_research_samples_file ${keep_research_samples_file}" : ""} \
        ${exclude_samples_file ? "--exclude_samples_file ${exclude_samples_file}" : ""}
    """

}
