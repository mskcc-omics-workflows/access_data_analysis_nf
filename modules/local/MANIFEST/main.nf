// manifest-cli's image lives in a private JFrog registry that Nextflow can't
// authenticate to on its own (pulling a Singularity image happens in
// Nextflow's own head process, before any task exists, so there's no process
// directive that can supply registry credentials for it -- see the MANIFEST
// module's git history if that's ever revisited). So the image is pre-seeded
// into Nextflow's own Singularity cache dir instead of pulled live -- the
// container reference below is normal, but the pull step is manual:
//   singularity pull \
//     /data1/core005/voyager_staging/ridgeback/voyager_pipeline_cache/mskcc.jfrog.io-omicswf-docker-prod-local-manifest_cli-0.1.0.img \
//     docker://mskcc.jfrog.io/omicswf-docker-prod-local/manifest_cli:0.1.0
// That filename isn't arbitrary -- it's Nextflow's own cache-key derivation
// (image ref with '/' and ':' replaced by '-', plus .img), so it has to match
// exactly or Nextflow treats it as a miss and attempts a live pull again.
// Re-run that (bumping the tag/filename) whenever manifest_cli is
// rebuilt/retagged -- nothing here will pick up a new image automatically.
process MANIFEST {
    tag "$request_id"
    label 'process_single'

    container 'mskcc.jfrog.io/omicswf-docker-prod-local/manifest_cli:0.1.0'

    input:
    val request_id
    path manifest_secrets_file

    publishDir "${params.outdir}/manifest", mode: params.publish_dir_mode

    output:
    path "${request_id}/manifest.tsv", emit: manifest_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    // manifest-cli's own credential cache/login is unused here: DATABRICKS_*
    // are exported straight into the environment, which config.py already
    // prefers over the ~/.manifest cache. --output pins the file inside the
    // task dir regardless of what MANIFEST_OUTPUT_DIR the secrets file sets.
    """
    set -a
    source ${manifest_secrets_file}
    set +a

    manifest generate ${request_id} --output .
    """
}
