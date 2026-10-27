/*
 * Uses the R library 'scuttle' to aggregate scRNA-seq to pseudobulk by sample and guide information.
 */

process pseudobulk {
    
    input:
    path data
    val celltype
    val facLoc
    val nOn
    val koEff
    val nGuide
    
    output:
    tuple path("celltype_${celltype}_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}.h5ad")
          val(celltype),
          val(facLoc), 
          val(nOn),
          val(koEff),
          val(nGuide)
    
    script:
    """
    Rscript ${workflow.launchDir}/R/pseudobulk.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${celltype}" \
        "${facLoc}" \
        "${nOn}" \
        "${koEff}" \
        "${nGuide}"
    """
}