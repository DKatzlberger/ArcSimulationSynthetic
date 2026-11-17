/*
 * Takes a Seurat object and prepares it for simulation permutation. It subsets to the defined celltype and knockout.
 */

params.script_id = "01"

process prepare {

    input:
    path data
    val celltype
    val knockout

    output:
    path("${params.script_id}-celltype_${celltype}_knockout_${knockout}.h5ad")

    script:
    """
    Rscript ${workflow.launchDir}/R/01-prepare.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${celltype}" \
        "${knockout}"
    """
}