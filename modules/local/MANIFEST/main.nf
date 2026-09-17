// manifest-cli's image lives in a private JFrog registry, so pulling it needs
// auth. Every profile in this repo runs under Singularity (voyager/juno/iris,
// and even the "docker" profile flips singularity.enabled = true), and
// Nextflow pulls Singularity images itself in its own head process -- before
// any task or work dir exists -- reading SINGULARITY_DOCKER_USERNAME /
// SINGULARITY_DOCKER_PASSWORD from the environment `nextflow run` was
// launched in. That's outside anything a process directive (beforeScript,
// path inputs) can reach, so there's nothing to wire in here: whoever
// launches the pipeline must `source` the secrets file (or otherwise export
// those two vars) before invoking `nextflow run`.
process MANIFEST {
    tag "$request_id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'mskcc.jfrog.io/omicswf-docker-prod-local/manifest_cli:0.1.0':
        'mskcc.jfrog.io/omicswf-docker-prod-local/manifest_cli:0.1.0' }"

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
