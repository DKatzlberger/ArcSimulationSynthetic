/*
 * Uses the R library 'limma' to do DGE analysis.
 */

process limma {

    input:
    path data
    val facLoc
    val nOn
    val koEff
    val nGuide
    val formula
    val ref
    
    output:
    tuple path("limma_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}_${formula}_fit.rds"), 
          path("limma_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/limma.R \
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