// Modules
include { splatter }   from "../modules/splatter.nf"
include { scdesign3_offtarget_mixed } from "../modules/scdesign3_offtarget_mixed.nf"
include { subset }     from "../modules/subset.nf"
include { pseudobulk } from "../modules/pseudobulk.nf"
include { limma }      from "../modules/limma.nf"
include { edger }      from "../modules/edger.nf"
include { deseq2 }     from "../modules/deseq2.nf"
include { dream }      from "../modules/dream.nf"
include { nebula }     from "../modules/nebula.nf"
include { sceptre }    from "../modules/sceptre.nf"
include { wilcox }     from "../modules/wilcox.nf"
include { scmageck }   from "../modules/scmageck.nf"

// Parameters
params {
    data: Path = "data/GSE213511/rna_integrated.rds"
    
    // splatter
    splatter_celltype = ["Ery"]
    splatter_facLoc   = [ 0.1, 0.5, 1.5, 1, 2 ]
    
    // scdesign3
    scdesign3_nOn   = [ 10, 50, 100, 200, 500 ]
    scdesign3_koEff = [ 25, 50, 75, 100 ]

    // subset
    subset_nGuide = [ 10, 20, 50, 100, 200, 500, 1000 ]
    
    // limma
    limma_formula = [ "~genotype" , "~genotype + sample" , "~guide" , "~guide + sample" ]
    limma_ref     = [ "NTC" , "NTC" ,  "NTC-1" , "NTC-1" ]
    
    // edger
    edger_formula = [ "~genotype" , "~genotype + sample" , "~guide" , "~guide + sample" ]
    edger_ref     = [ "NTC" , "NTC" ,  "NTC-1" , "NTC-1" ]

    // deseq2
    deseq2_formula = [ "~genotype" , "~genotype + sample" , "~guide" , "~guide + sample" ]
    deseq2_ref     = [ "NTC" , "NTC" ,  "NTC-1" , "NTC-1" ]
    
    // dream
    dream_formula = [ "~genotype" , "~genotype + sample" , "~genotype + (1|sample)" , "~genotype + sample + (1|guide)", "~genotype + (1|sample) + (1|guide)" , "~guide" , "~guide + sample" , "~guide + (1|sample)" ]
    dream_ref     = [ "NTC" , "NTC" , "NTC" , "NTC", "NTC" , "NTC-1" , "NTC-1" , "NTC-1" ]
    
    // nebula
    nebula_formula = [ "~genotype" , "~genotype + sample" , "~genotype + (1|sample)" , "~genotype + sample + (1|guide)", "~guide" , "~guide + sample" , "~guide + (1|sample)" ]
    nebula_ref     = [ "NTC" , "NTC" , "NTC" , "NTC", "NTC-1" , "NTC-1" , "NTC-1" ]
    
    // sceptre
    sceptre_formula = [ "~genotype" , "~genotype + sample" ]
    sceptre_ref     = [ "NTC" , "NTC" ]
    
    // wilcox
    wilcox_formula = [ "~genotype" ]
    wilcox_ref     = [ "NTC" ]

    // scmageck
    scmageck_formula = [ "~genotype" ]
    scmageck_ref     = [ "NTC" ]
}

workflow {
    main:

    // Simualtions ============================================================
    // R: splatter ------------------------------------------------------------
    splatter_args = channel.value(params.data)
        .combine(channel.fromList(params.splatter_celltype))
        .combine(channel.fromList(params.splatter_facLoc))

    splatter(
        splatter_args.map { arg -> arg[0] }, // data
        splatter_args.map { arg -> arg[1] }, // celltype
        splatter_args.map { arg -> arg[2] }  // facLoc
    )

    // R: scDesign3 -----------------------------------------------------------
    scdesign3_args = splatter.out
        .combine(params.scdesign3_nOn)
        .combine(params.scdesign3_koEff)

    scdesign3_offtarget_mixed(
        scdesign3_args.map { arg -> arg[0] }, // data 
        scdesign3_args.map { arg -> arg[1] }, // celltype 
        scdesign3_args.map { arg -> arg[2] }, // facLoc 
        scdesign3_args.map { arg -> arg[3] }, // nOn
        scdesign3_args.map { arg -> arg[4] }  // koEff
    )

    // R: subset --------------------------------------------------------------
    subset_args = scdesign3_offtarget_mixed.out
        .combine(channel.fromList(params.subset_nGuide))

    subset(
        subset_args.map { arg -> arg[0] }, // data 
        subset_args.map { arg -> arg[2] }, // celltype 
        subset_args.map { arg -> arg[3] }, // facLoc 
        subset_args.map { arg -> arg[4] }, // nOn
        subset_args.map { arg -> arg[5] }, // koEff
        subset_args.map { arg -> arg[6] }  // nGuide
    )

    // R: pseudobulk ----------------------------------------------------------
    pseudobulk_args = subset.out

    pseudobulk(
        pseudobulk_args.map { arg -> arg[0] }, // data
        pseudobulk_args.map { arg -> arg[1] }, // celltype
        pseudobulk_args.map { arg -> arg[2] }, // facLoc
        pseudobulk_args.map { arg -> arg[3] }, // nOn
        pseudobulk_args.map { arg -> arg[4] }, // koEff
        pseudobulk_args.map { arg -> arg[5] }  // nGuide
    )

    // DGE pseudobulk =========================================================
    // R: limma ---------------------------------------------------------------
    limma_formula = channel.fromList(params.limma_formula)
    limma_ref     = channel.fromList(params.limma_ref)
    limma_pairs   = limma_formula.merge(limma_ref)

    limma_args = pseudobulk.out.combine(limma_pairs)

    limma(
        limma_args.map { arg -> arg[0] }, // data
        limma_args.map { arg -> arg[1] }, // celltype
        limma_args.map { arg -> arg[2] }, // facLoc
        limma_args.map { arg -> arg[3] }, // nOn
        limma_args.map { arg -> arg[4] }, // koEff
        limma_args.map { arg -> arg[5] }, // nGuide
        limma_args.map { arg -> arg[6] }, // formula
        limma_args.map { arg -> arg[7] }  // ref
    )

    // R: edger ---------------------------------------------------------------
    edger_formula = channel.fromList(params.edger_formula)
    edger_ref     = channel.fromList(params.edger_ref)
    edger_pairs   = edger_formula.merge(edger_ref)

    edger_args = pseudobulk.out.combine(edger_pairs)

    edger(
        edger_args.map { arg -> arg[0] }, // data
        edger_args.map { arg -> arg[1] }, // celltype
        edger_args.map { arg -> arg[2] }, // facLoc
        edger_args.map { arg -> arg[3] }, // nOn
        edger_args.map { arg -> arg[4] }, // koEff
        edger_args.map { arg -> arg[5] }, // nGuide
        edger_args.map { arg -> arg[6] }, // formula
        edger_args.map { arg -> arg[7] }  // ref
    )

    // R: deseq2 --------------------------------------------------------------
    deseq2_formula = channel.fromList(params.deseq2_formula)
    deseq2_ref     = channel.fromList(params.deseq2_ref)
    deseq2_pairs   = deseq2_formula.merge(deseq2_ref)

    deseq2_args = pseudobulk.out.combine(deseq2_pairs)

    deseq2(
        deseq2_args.map { arg -> arg[0] }, // data
        deseq2_args.map { arg -> arg[1] }, // celltype
        deseq2_args.map { arg -> arg[2] }, // facLoc
        deseq2_args.map { arg -> arg[3] }, // nOn
        deseq2_args.map { arg -> arg[4] }, // koEff
        deseq2_args.map { arg -> arg[5] }, // nGuide
        deseq2_args.map { arg -> arg[6] }, // formula
        deseq2_args.map { arg -> arg[7] }  // ref
    )

    // R: dream ---------------------------------------------------------------
    dream_formula = channel.fromList(params.dream_formula)
    dream_ref     = channel.fromList(params.dream_ref)
    dream_pairs   = dream_formula.merge(dream_ref)

    dream_args = pseudobulk.out.combine(dream_pairs)

    dream(
        dream_args.map { arg -> arg[0] }, // data
        dream_args.map { arg -> arg[1] }, // celltype
        dream_args.map { arg -> arg[2] }, // facLoc
        dream_args.map { arg -> arg[3] }, // nOn
        dream_args.map { arg -> arg[4] }, // koEff
        dream_args.map { arg -> arg[5] }, // nGuide
        dream_args.map { arg -> arg[6] }, // formula
        dream_args.map { arg -> arg[7] }  // ref
    )

    // DGE single cell ========================================================
    // R: nebula --------------------------------------------------------------
    nebula_formula = channel.fromList(params.nebula_formula)
    nebula_ref     = channel.fromList(params.nebula_ref)
    nebula_pairs   = nebula_formula.merge(nebula_ref)

    nebula_args = subset.out.combine(nebula_pairs)

    nebula(
        nebula_args.map { arg -> arg[0] }, // data
        nebula_args.map { arg -> arg[1] }, // celltype
        nebula_args.map { arg -> arg[2] }, // facLoc
        nebula_args.map { arg -> arg[3] }, // nOn
        nebula_args.map { arg -> arg[4] }, // koEff
        nebula_args.map { arg -> arg[5] }, // nGuide
        nebula_args.map { arg -> arg[6] }, // formula
        nebula_args.map { arg -> arg[7] }  // ref
    )

    // R: sceptre -------------------------------------------------------------
    sceptre_formula = channel.fromList(params.sceptre_formula)
    sceptre_ref     = channel.fromList(params.sceptre_ref)
    sceptre_pairs   = sceptre_formula.merge(sceptre_ref)

    sceptre_args = subset.out.combine(sceptre_pairs)

    sceptre(
        sceptre_args.map { arg -> arg[0] }, // data
        sceptre_args.map { arg -> arg[1] }, // celltype
        sceptre_args.map { arg -> arg[2] }, // facLoc
        sceptre_args.map { arg -> arg[3] }, // nOn
        sceptre_args.map { arg -> arg[4] }, // koEff
        sceptre_args.map { arg -> arg[5] }, // nGuide
        sceptre_args.map { arg -> arg[6] }, // formula
        sceptre_args.map { arg -> arg[7] }  // ref
    )

    // R: wilcox --------------------------------------------------------------
    wilcox_formula = channel.fromList(params.wilcox_formula)
    wilcox_ref     = channel.fromList(params.wilcox_ref)
    wilcox_pairs   = wilcox_formula.merge(wilcox_ref)

    wilcox_args = subset.out.combine(wilcox_pairs)

    wilcox(
        wilcox_args.map { arg -> arg[0] }, // data
        wilcox_args.map { arg -> arg[1] }, // celltype
        wilcox_args.map { arg -> arg[2] }, // facLoc
        wilcox_args.map { arg -> arg[3] }, // nOn
        wilcox_args.map { arg -> arg[4] }, // koEff
        wilcox_args.map { arg -> arg[5] }, // nGuide
        wilcox_args.map { arg -> arg[6] }, // formula
        wilcox_args.map { arg -> arg[7] }  // ref
    )

    // R: scmageck ------------------------------------------------------------
    scmageck_formula = channel.fromList(params.scmageck_formula)
    scmageck_ref     = channel.fromList(params.scmageck_ref)
    scmageck_pairs   = scmageck_formula.merge(scmageck_ref)

    scmageck_args = subset.out.combine(scmageck_pairs)

    scmageck(
        scmageck_args.map { arg -> arg[0] }, // data
        scmageck_args.map { arg -> arg[1] }, // celltype
        scmageck_args.map { arg -> arg[2] }, // facLoc
        scmageck_args.map { arg -> arg[3] }, // nOn
        scmageck_args.map { arg -> arg[4] }, // koEff
        scmageck_args.map { arg -> arg[5] }, // nGuide
        scmageck_args.map { arg -> arg[6] }, // formula
        scmageck_args.map { arg -> arg[7] }  // ref
    )

    publish:
    // Sim. res
    res_splatter   = splatter.out
    res_scdesign3_offtarget_mixed  = scdesign3_offtarget_mixed.out
    res_subset     = subset.out
    res_pseudobulk = pseudobulk.out

    // DGE res
    res_limma      = limma.out
    res_edger      = edger.out
    res_deseq2     = deseq2.out
    res_dream      = dream.out
    res_nebula     = nebula.out
    res_sceptre    = sceptre.out
    res_wilcox     = wilcox.out
    res_scmageck   = scmageck.out
}

output {
    res_splatter {
        path "nextflow/ArcSimSyn_offtarget_mixed/splatter"
        mode "copy"
    }
    scdesign3_offtarget_mixed {
        path "nextflow/ArcSimSyn_offtarget_mixed/scdesign3"
        mode "copy"
    }
    res_subset {
        path "nextflow/ArcSimSyn_offtarget_mixed/subset"
        mode "copy"
    }
    res_pseudobulk {
        path "nextflow/ArcSimSyn_offtarget_mixed/pseudobulk"
        mode "copy"
    }

    res_limma {
        path "nextflow/ArcSimSyn_offtarget_mixed/limma"
        mode "copy"
    }
    res_edger {
        path "nextflow/ArcSimSyn_offtarget_mixed/edger"
        mode "copy"
    }
    res_deseq2 {
        path "nextflow/ArcSimSyn_offtarget_mixed/deseq2"
        mode "copy"
    }
    res_dream {
        path "nextflow/ArcSimSyn_offtarget_mixed/dream"
        mode "copy"
    }
    res_nebula {
        path "nextflow/ArcSimSyn_offtarget_mixed/nebula"
        mode "copy"
    }
    res_sceptre {
        path "nextflow/ArcSimSyn_offtarget_mixed/sceptre"
        mode "copy"
    }
    res_wilcox {
        path "nextflow/ArcSimSyn_offtarget_mixed/wilcox"
        mode "copy"
    }
    res_scmageck {
        path "nextflow/ArcSimSyn_offtarget_mixed/scmageck"
        mode "copy"
    }
}