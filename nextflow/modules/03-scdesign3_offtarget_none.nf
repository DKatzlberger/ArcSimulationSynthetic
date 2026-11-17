params.script_id = "03"

process scdesign3_offtarget_none {
    
    input:
    path data
    val fracOff
    val koEff

    output:
    tuple path("${params.script_id}-${data.baseName}_fracOff_${fracOff}_koEff_${koEff}.h5ad"),
          path("${params.script_id}-${data.baseName}_fracOff_${fracOff}_koEff_${koEff}_A.svg"),
          path("${params.script_id}-${data.baseName}_fracOff_${fracOff}_koEff_${koEff}_B.svg"),
          path("${params.script_id}-${data.baseName}_fracOff_${fracOff}_koEff_${koEff}_C.svg")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/03-scdesign3_offtarget_none.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${fracOff}" \
        "${koEff}"
    """
}