
/*
 * Use the R library 'splatter' to simulate single cell data with two groups and 3 batches.
 */

process splatter {
    
    input:
    path data
    val facLoc
    
    output:
    tuple val(facLoc),
          path("splatter_facLoc_${facLoc}.h5ad"), 
          path("splatter_facLoc_${facLoc}_bx.svg"), 
          path("splatter_facLoc_${facLoc}_ht.svg"), 
          path("splatter_facLoc_${facLoc}_um.svg")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/splatter.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${facLoc}"
    """
}