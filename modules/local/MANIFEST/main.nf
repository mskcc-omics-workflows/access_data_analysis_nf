// manifest-cli's image lives in a private JFrog registry, so the docker daemon
// needs to be logged in before it can pull it. beforeScript runs on the host,
// before the container is launched, so it's the one place that can do that login.
// We skip it for singularity since there's no docker daemon to log in.
def jfrog_login_cmd = {
    if (workflow.containerEngine != 'docker') {
        return ''
    }
    "set -a; source ${params.manifest_secrets_file}; set +a; echo \"\$JFROG_TOKEN\" | docker login mskcc.jfrog.io -u \"\$JFROG_USER\" --password-stdin"
}

process MANIFEST {
    tag "$request_id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'mskcc.jfrog.io/omicswf-docker-prod-local/manifest_cli:0.1.0':
        'mskcc.jfrog.io/omicswf-docker-prod-local/manifest_cli:0.1.0' }"

    beforeScript jfrog_login_cmd()

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
