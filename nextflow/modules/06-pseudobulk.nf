/*
 * Uses the R library 'scuttle' to aggregate scRNA-seq to pseudobulk by sample and guide information.
 */

params.script_id = "06"

process pseudobulk {
    
    input:
    path data
    val groups
    
    output:
    path("${params.script_id}-${data.baseName}.h5ad")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/06-pseudobulk.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${groups}"
    """
}