/*
 * Use the R library 'scdesign3' to recreate scRNA-seq dataset. 
 * Two unique gene sets are conditioned on gene expression from a knockout.
 * One gene set is the 'real' downstream program of the knockout.
 * The second gene set are off-target downstream genes from an off-target guide.
 * The gene sets overlap but the off-target guides has additional (more) downstream genes compared to on-target.
 */

process scdesign3_offtarget_mixed {
    
    input:
    path data
    val celltype
    val facLoc
    val nOn
    val koEff

    output:
    tuple path("celltype_${celltype}_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}.h5ad"),
          val(celltype)
          val(facLoc),
          val(nOn),
          val(koEff), 
          path("celltype_${celltype}_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_bx.svg"), 
          path("celltype_${celltype}_facLoc_${facLoc}_nOn_${nOn}_koEff_${koEff}_ht.svg")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/scdesign3_offtarget_mixed.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${celltype}" \
        "${facLoc}" \
        "${nOn}" \
        "${koEff}"
    """
}