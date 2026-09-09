/*
 * Uses the R library 'sceptre' to do DGE analysis.
 */

process sceptre {

    input:
    path data
    val facLoc
    val nOn
    val koEff
    val nGuide
    val formula
    val ref
    
    output:
    tuple path("sceptre_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}_${formula}_fit.rds"), 
          path("sceptre_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/sceptre.R \
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