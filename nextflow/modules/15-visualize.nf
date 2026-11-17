params.script_id = "15"

process visualize {

    input:
    path dge_files

    output:
    path "panel_*", emit: panels


    script:
    """
    Rscript ${workflow.launchDir}/R/15-visualize.R \
        "${workflow.launchDir}" \
        ${dge_files}
    """
}