process INFER_SAMPLES {
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path id_mapping_file
    path include_samples_file
    path exclude_samples_file
    path dmp_access_key_path, name: "access_key.txt"
    path dmp_impact_key_path, name: "dmp_key.txt"
    val research_access_dir_template
    val access_sample_regex_pattern
    val impact_sample_regex_pattern
        

    publishDir '${outdir}/intermediary/patient_JSONs', mode: 'copy'

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
        --dmp_access_key_path $dmp_access_key_path \\
        --dmp_impact_key_path $dmp_impact_key_path \\
        --research_access_dir_template $research_access_dir_template \\
        --access_sample_regex_pattern '$access_sample_regex_pattern' \\
        --impact_sample_regex_pattern '$impact_sample_regex_pattern' \\
    """

}