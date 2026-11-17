
params.script_id = "04"

process compare {

    input:
    path data
    path real_data
    val celltype

    output:
    path("${params.script_id}-${data.baseName}_A.svg")


    script:
    """
    Rscript ${workflow.launchDir}/R/04-compare.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${real_data}" \
        "${celltype}"
    """
}