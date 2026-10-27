/*
 * Uses the R library 'limma' to do DGE analysis.
 */

process limma {

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
    Rscript ${workflow.launchDir}/R/limma.R \
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