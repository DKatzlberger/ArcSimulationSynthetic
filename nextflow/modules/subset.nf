/*
 * Subsamples the KO cells to reduce the number of cells with the knockout.
 * Keeps all of the NTC cells.
 */

process subset {
    
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
          path("subset_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}.h5ad")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/subset.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${facLoc}" \
        "${nOn}" \
        "${koEff}" \
        "${nGuide}"
    """
}