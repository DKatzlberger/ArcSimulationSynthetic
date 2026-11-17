/*
 * Subsamples the KO cells to reduce the number of cells with the knockout.
 * Keeps all of the NTC cells.
 */

params.script_id = "05"

process subset {
    
    input:
    path data
    val nGuide

    output:
    path("${params.script_id}-${data.baseName}_nGuide_${nGuide}.h5ad")

    script:
    """
    Rscript ${workflow.launchDir}/R/05-subset.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${nGuide}"
    """
}