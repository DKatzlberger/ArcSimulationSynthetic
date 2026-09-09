/*
 * Uses the R library 'Seurat' to do DGE analysis.
 */

process wilcox {

    input:
    path data
    val facLoc
    val nOn
    val koEff
    val nGuide
    val formula
    val ref
    
    output:
    tuple path("wilcox_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}_${formula}_fit.rds"), 
          path("wilcox_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/wilcox.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${facLoc}" \
        "${nOn}" \
        "${koEff}" \
        "${nGuide}" \
        "${formula}" \
        "${ref}"
    """
}