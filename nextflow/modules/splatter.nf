/*
 * Use the R library 'splatter' to simulate single cell data with two groups and 3 batches.
 */

process splatter {
    
    input:
    path data
    val celltype
    val facLoc
    
    output:
    tuple path("celltype_${celltype}_facLoc_${facLoc}.h5ad"),
          val(celltype),
          val(facLoc), 
          path("celltype_${celltype}_facLoc_${facLoc}_bx.svg"), 
          path("celltype_${celltype}_facLoc_${facLoc}_ht.svg")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/splatter.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${celltype}" \
        "${facLoc}"
    """
}