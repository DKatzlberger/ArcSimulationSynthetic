/*
 * Uses the R library 'DESeq2' to do DGE analysis.
 */

process deseq2 {

    input:
    path data
    val facLoc
    val nOn
    val koEff
    val nGuide
    val formula
    val ref
    
    output:
    tuple path("deseq2_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}_${formula}_fit.rds"), 
          path("deseq2_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/deseq2.R \
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