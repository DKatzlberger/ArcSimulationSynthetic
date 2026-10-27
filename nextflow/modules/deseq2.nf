/*
 * Uses the R library 'DESeq2' to do DGE analysis.
 */

process deseq2 {

    input:
    path data
    val celltype
    val facLoc
    val nOn
    val koEff
    val nGuide
    val formula
    val ref
    
    output:
    path("celltype_${celltype}_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/deseq2.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${celltype}" \
        "${facLoc}" \
        "${nOn}" \
        "${koEff}" \
        "${nGuide}" \
        "${formula}" \
        "${ref}"
    """
}