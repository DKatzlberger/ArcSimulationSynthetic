/*
 * Use the R library 'scdesign3' to recreate scRNA-seq dataset. 
 * Two unique gene sets are conditioned on gene expression from a knockout.
 * One gene set is the 'real' downstream program of the knockout.
 * The second gene set are off-target downstream genes from an off-target guide.
 * The gene sets don't overlap.
 */

process scdesign3 {
    
    input:
    path data
    val facLoc
    val nOn
    val koEff

    output:
    tuple val(facLoc), 
          val(nOn),
          val(koEff),
          path("scdesign3_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}.h5ad"), 
          path("scdesign3_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_bx.svg"), 
          path("scdesign3_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_ht.svg"), 
          path("scdesign3_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_um.svg")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/scdesign3_offtarget_unique.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${facLoc}" \
        "${nOn}" \
        "${koEff}"
    """
}