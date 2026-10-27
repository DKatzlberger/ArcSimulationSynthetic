/*
 * Subsamples the KO cells to reduce the number of cells with the knockout.
 * Keeps all of the NTC cells.
 */

process subset {
    
    input:
    path data
    val celltype
    val facLoc
    val nOn
    val koEff
    val nGuide

    output:
    tuple path("celltype_${celltype}_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_nGuide_${nGuide}.h5ad"),
          val(celltype),
          val(facLoc),
          val(nOn),
          val(koEff),
          val(nGuide)
    
    script:
    """
    Rscript ${workflow.launchDir}/R/subset.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${celltype}" \
        "${facLoc}" \
        "${nOn}" \
        "${koEff}" \
        "${nGuide}"
    """
}