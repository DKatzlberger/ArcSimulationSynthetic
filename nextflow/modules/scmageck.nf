/*
 * Uses the R library 'scMAGeCK' to do DGE analysis.
 */

process scmageck {

    input:
    path data
    val facLoc
    val nOn
    val koEff
    val nGuide
    val formula
    val ref
    
    output:
    tuple path("scmageck_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}_${formula}_fit.rds"), 
          path("scmageck_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/scmageck.R \
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