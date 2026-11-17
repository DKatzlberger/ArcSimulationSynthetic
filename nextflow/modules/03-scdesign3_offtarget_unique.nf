/*
 * Use the R library 'scdesign3' to recreate scRNA-seq dataset. 
 * Two unique gene sets are conditioned on gene expression from a knockout.
 * One gene set is the 'real' downstream program of the knockout.
 * The second gene set are off-target downstream genes from an off-target guide.
 * The gene sets don't overlap.
 */

params.script_id = "03"

process scdesign3_offtarget_unique {
    
    input:
    path data
    val nOn
    val koEff

    output:
    tuple path("${params.script_id}-${data.baseName}_nOn_${nOn}_koEff_${koEff}.h5ad"),
          path("${params.script_id}-${data.baseName}_nOn_${nOn}_koEff_${koEff}_ht.svg"),
          path("${params.script_id}-${data.baseName}_nOn_${nOn}_koEff_${koEff}_bx.svg")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/03-scdesign3_offtarget_unique.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${nOn}" \
        "${koEff}"
    """
}