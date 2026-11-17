/*
 * Use the R library 'splatter' to simulate single cell data with two groups and 3 batches.
 */

params.script_id = "02"

process splatter {
    
    input:
    path data
    val facLoc
    
    output:
    tuple path("${params.script_id}-${data.baseName}_facLoc_${facLoc}.h5ad"),
          path("${params.script_id}-${data.baseName}_facLoc_${facLoc}_A.svg"),
          path("${params.script_id}-${data.baseName}_facLoc_${facLoc}_B.svg"),
          path("${params.script_id}-${data.baseName}_facLoc_${facLoc}_C.svg")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/02-splatter.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${facLoc}"
    """
}