/*
 * Uses the R library 'scuttle' to aggregate scRNA-seq to pseudobulk by sample and guide information.
 */

process pseudobulk {
    
    input:
    path data
    val facLoc
    val nOn
    val koEff
    val nGuide
    
    output:
    tuple val(facLoc), 
          val(nOn),
          val(koEff),
          val(nGuide),
          path("pseudobulk_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}.h5ad")

    
    script:
    """
    Rscript ${workflow.launchDir}/R/pseudobulk.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${facLoc}" \
        "${nOn}" \
        "${koEff}" \
        "${nGuide}"
    """
}