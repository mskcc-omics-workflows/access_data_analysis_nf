/*
* Step 1: Aggregate allele counts across all samples
* - Get allele counts from genotyped MAFs and merge them with union calls
* - Applies coverage thresholds (ACCESS and IMPACT specific)
* - Calculates VAF and determines call status for each variant
* - call status is either called, genotyped or low_coverage
* - Creates a comprehensive table of variants across all samples
*/
process SNV_INDEL_AGGREGATE_ALLELE_COUNTS {
    label 'process_single'
    tag "${patient_id}"

    input:
        tuple path(patient_json), val(patient_id), path(genotyping_output), path(all_small_calls_maf)
        val access_min_cov
        val impact_min_cov

    publishDir "${params.outdir}/intermediary/small_variants/${patient_id}", mode: 'copy', pattern: '*SNV-INDEL.allele_counts.csv'

    output:
        tuple val(patient_id), path("*SNV-INDEL.allele_counts.csv"), emit: snv_indel_results

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    python3 ${workflow.projectDir}/bin/aggregate_snv_indel_allele_counts.py \\
        --patient_json $patient_json \\
        --genotyped_mafs $genotyping_output \\
        --union_calls_maf $all_small_calls_maf \\
        --output ${patient_id}-SNV-INDEL.allele_counts.csv \\
        --access_min_cov ${access_min_cov} \\
        --impact_min_cov ${impact_min_cov}
    """
}

/*
* Step 2: Annotate variants with Hotspot and CH status
* - Compares variants against known hotspot mutations
* - Identifies potential clonal hematopoiesis (CH) variants
* - Adds two new columns: 'Hotspot' and 'CH' (yes/"" )
* - Maintains all original variant information
*/
process SNV_INDEL_ANNOTATE_HOTSPOT_CH {
    label 'process_single'
    tag "${patient_id}"

    input:
        tuple val(patient_id), path(snv_indel_csv)
        val hotspot_list
        val ch_list

    publishDir "${params.outdir}/intermediary/small_variants/${patient_id}", mode: 'copy', pattern: "*hotspot_ch.csv"
    output:
        tuple val(patient_id), path("*hotspot_ch.csv"), emit: hotspot_ch_annotated_snv_indel

    when:
        task.ext.when == null || task.ext.when

    script:
    def basename = snv_indel_csv.baseName
    """
    python3 ${workflow.projectDir}/bin/annotate_snv_indel_hotspot_ch.py \\
        --variant_input ${snv_indel_csv} \\
        --hotspot_list ${hotspot_list} \\
        --ch_list ${ch_list} \\
        --output ${basename}.hotspot_ch.csv
    """
}

/*
* Step 3: Add filter annotations
* - Adds a 'filter' column to indicate variant status
* - Marks variants from excluded genes
* - Marks variants with excluded classifications
* - Marks variants with low coverage in ACCESS samples
* - Combines multiple filter reasons with semicolons
* - Produces both intermediate output with all rows
*   and final filtered output with only passing variants
*/
process SNV_INDEL_ADD_FILTER_COL {
    label 'process_single'
    tag "${patient_id}"

    input:
        tuple val(patient_id), path(snv_indel_csv)
        val exclude_genes
        val exclude_classifications

    publishDir "${params.outdir}/intermediary/small_variants/${patient_id}", mode: 'copy', pattern: "*filter.csv"
    publishDir "${params.outdir}/final/${patient_id}", mode: 'copy', pattern: "*snv_indel.filtered.csv"
    output:
        tuple val(patient_id), path("*filter.csv"), emit: filtered_snv_indel
        tuple val(patient_id), path("*snv_indel.filtered.csv"), emit: final_filtered_snv_indel

    when:
        task.ext.when == null || task.ext.when

    script:
    def basename = snv_indel_csv.baseName
    """
    python3 ${workflow.projectDir}/bin/filter_snv_indel_calls.py \\
        --variant_input ${snv_indel_csv} \\
        --exclude_genes ${exclude_genes} \\
        --exclude_classifications ${exclude_classifications} \\
        --output ${basename}.filter.csv \\
        --output_final ${patient_id}.snv_indel.filtered.csv
    """
}

/*
* Step 4: Calculate FACETS-adjusted VAF
* - Processes variants with FACETS copy number data
* - Calculates adjusted VAF considering copy number changes
* - Only adjusts VAF for clonal variants
* - Handles multiple FACETS IMPACT samples per patient
* - Selects best FACETS IMPACT sample based on most overlap with ACCESS variants
* - Produces files with adjusted VAFs using all FACETS IMPACT samples,
*   and best FACETS IMPACT sample
*/
process SNV_INDEL_ADD_FACETS_ADJUSTED_VAF {
    label 'process_single'
    tag "${patient_id}"

    input:
        tuple val(patient_id), path(facets_fit), path(snv_indel_csv)

    publishDir "${params.outdir}/intermediary/small_variants/${patient_id}", mode: 'copy', pattern: '*adj_vaf_all_impact.csv'
    publishDir "${params.outdir}/final/${patient_id}", mode: 'copy', pattern: '*snv_indel.filtered.adj_vaf.csv'
    output:
        tuple val(patient_id), path("*adj_vaf*.csv"), emit: adjusted_vaf_results

    when:
        task.ext.when == null || task.ext.when

    script:
    def basename = snv_indel_csv.baseName
    """
    python3 ${workflow.projectDir}/bin/add_facets_adjusted_vaf.py \\
        --variant_csv ${snv_indel_csv} \\
        --facets_file_list ${facets_fit} \\
        --output_all_facets_samples ${basename}.adj_vaf_all_impact.csv \\
        --output_best_facets_sample ${patient_id}.snv_indel.filtered.adj_vaf.csv
    """
}


