process INFER_SAMPLES {
    label 'process_single'

    input:
    path id_mapping_file
    path include_samples_file
    path exclude_samples_file
    path clinical_access_key_file, name: "access_key.txt"
    path clinical_impact_key_file, name: "dmp_key.txt"
    val research_access_bam_dir
    val clinical_access_sample_regex_pattern
    val clinical_impact_sample_regex_pattern        

    publishDir "${params.outdir}/intermediary/patient_JSONs", mode: 'copy'

    output:
    path "*.json", emit: all_samples_json


    when:
    task.ext.when == null || task.ext.when

    script:

    """
    python3 ../../../bin/infer_samples.py \\
        --id_mapping_file $id_mapping_file \\
        --include_samples_file $include_samples_file \\
        --exclude_samples_file $exclude_samples_file \\
        --clinical_access_key_file $clinical_access_key_file \\
        --clinical_impact_key_file $clinical_impact_key_file \\
        --research_access_bam_dir $research_access_bam_dir \\
        --clinical_access_sample_regex_pattern '$clinical_access_sample_regex_pattern' \\
        --clinical_impact_sample_regex_pattern '$clinical_impact_sample_regex_pattern' \\
    """

}