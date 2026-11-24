// Modules
include { prepare }                   from "../modules/01-prepare.nf"
include { splatter }                  from "../modules/02-splatter.nf"
include { scdesign3_offtarget_unique } from "../modules/03-scdesign3_offtarget_unique.nf"
include { compare }                   from "../modules/04-compare.nf"
include { subset }                    from "../modules/05-subset.nf"
include { pseudobulk }                from "../modules/06-pseudobulk.nf"
include { limma }                     from "../modules/07-limma.nf"
include { edger }                     from "../modules/08-edger.nf"
include { deseq2 }                    from "../modules/09-deseq2.nf"
include { dream }                     from "../modules/10-dream.nf"
include { nebula }                    from "../modules/11-nebula.nf"
include { sceptre }                   from "../modules/12-sceptre.nf"
include { wilcox }                    from "../modules/13-wilcox.nf"
include { scmageck }                  from "../modules/14-scmageck.nf"
include { visualize }                 from "../modules/15-visualize.nf"


// Parameters
params {
    data: Path = "data_generated/GSE213511/rna_integrated.rds"

    // Prepare
    prepare_celltype = ["Gran. P"]
    prepare_knockout = ["NULL"]
    
    // splatter
    splatter_facLoc = [ 0.5, 1, 1.5, 2 ] // 0.5, 1, 1.5, 2 
    
    // scdesign3
    scdesign3_fracOff = [ 10, 20, 30, 40, 50 ]  // 10, 20, 30, 40, 50
    scdesign3_koEff   = [ 100 ]

    // subset
    subset_nGuide = [ 10, 20, 50, 100, 150, 200, 250 ] // 10, 20, 50, 100, 150, 200, 250

    // pseudobulk
    pseudobulk_grps = [ "sample,guide" ]
    
    // limma
    limma_formula = [ "~genotype" , "~genotype + sample" ]
    limma_ref     = [ "NTC" , "NTC" ]
    
    // edger
    edger_formula = [ "~genotype" , "~genotype + sample" ]
    edger_ref     = [ "NTC" , "NTC" ]

    // deseq2
    deseq2_formula = [ "~genotype" , "~genotype + sample" ]
    deseq2_ref     = [ "NTC" , "NTC" ]
    
    // dream
    dream_formula = [ "~genotype" , "~genotype + sample" , "~genotype + (1|sample)" , "~genotype + sample + (1|guide)", "~genotype + (1|sample) + (1|guide)" ]
    dream_ref     = [ "NTC" , "NTC" , "NTC" , "NTC", "NTC" ]
    
    // nebula
    nebula_formula = [ "~genotype" , "~genotype + sample" , "~genotype + (1|sample)" , "~genotype + sample + (1|guide)" ]
    nebula_ref     = [ "NTC" , "NTC" , "NTC" , "NTC" ]
    
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
    // R: 01-prepare ----------------------------------------------------------
    prepare_args = channel.value(params.data)
        .combine(channel.fromList(params.prepare_celltype))
        .combine(channel.fromList(params.prepare_knockout))

    prepare(
        prepare_args.map { arg -> arg[0] }, // data
        prepare_args.map { arg -> arg[1] }, // celltype
        prepare_args.map { arg -> arg[2] }  // knockout
    )

    // R: 02-splatter ---------------------------------------------------------
    splatter_args = prepare.out
        .combine(channel.fromList(params.splatter_facLoc))

    splatter(
        splatter_args.map { arg -> arg[0] }, // data
        splatter_args.map { arg -> arg[1] }  // facLoc
    )

    // R: 03-scdesign3 --------------------------------------------------------
    scdesign3_offtarget_unique_args = splatter.out
        .map { h5ad, _A, _B, _C -> h5ad }
        .combine(params.scdesign3_fracOff)
        .combine(params.scdesign3_koEff)

    scdesign3_offtarget_unique(
        scdesign3_offtarget_unique_args.map { arg -> arg[0] }, // data 
        scdesign3_offtarget_unique_args.map { arg -> arg[1] }, // nOn
        scdesign3_offtarget_unique_args.map { arg -> arg[2] }  // koEff
    )

    // R: 04-compare ----------------------------------------------------------
    compare_args = scdesign3_offtarget_unique.out
        .map { h5ad, _A, _B, _C -> h5ad }
        .combine(channel.value(params.data))
        .combine(channel.fromList(params.prepare_celltype))

    compare(
        compare_args.map { arg -> arg[0] }, // data
        compare_args.map { arg -> arg[1] }, // real_data
        compare_args.map { arg -> arg[2] }  // celltype
    )

    // R: 05-subset -----------------------------------------------------------
    subset_args = scdesign3_offtarget_unique.out
        .map { h5ad, _A, _B, _C -> h5ad }
        .combine(channel.fromList(params.subset_nGuide))

    subset(
        subset_args.map { arg -> arg[0] }, // data 
        subset_args.map { arg -> arg[1] }  // nGuide
    )

    // R: 06-pseudobulk -------------------------------------------------------
    pseudobulk_args = subset.out.combine(channel.fromList(params.pseudobulk_grps))

    pseudobulk(
        pseudobulk_args.map { arg -> arg[0] }, // data
        pseudobulk_args.map { arg -> arg[1] }  // groups
    )

    // DGE pseudobulk =========================================================
    // R: 07-limma ------------------------------------------------------------
    limma_formula = channel.fromList(params.limma_formula)
    limma_ref     = channel.fromList(params.limma_ref)
    limma_pairs   = limma_formula.merge(limma_ref)

    limma_args = pseudobulk.out.combine(limma_pairs)

    limma(
        limma_args.map { arg -> arg[0] }, // data
        limma_args.map { arg -> arg[1] }, // formula
        limma_args.map { arg -> arg[2] }  // ref
    )

    // R: 08-edger ------------------------------------------------------------
    edger_formula = channel.fromList(params.edger_formula)
    edger_ref     = channel.fromList(params.edger_ref)
    edger_pairs   = edger_formula.merge(edger_ref)

    edger_args = pseudobulk.out.combine(edger_pairs)

    edger(
        edger_args.map { arg -> arg[0] }, // data
        edger_args.map { arg -> arg[1] }, // formula
        edger_args.map { arg -> arg[2] }  // ref
    )

    // R: 09-deseq2 -----------------------------------------------------------
    deseq2_formula = channel.fromList(params.deseq2_formula)
    deseq2_ref     = channel.fromList(params.deseq2_ref)
    deseq2_pairs   = deseq2_formula.merge(deseq2_ref)

    deseq2_args = pseudobulk.out.combine(deseq2_pairs)

    deseq2(
        deseq2_args.map { arg -> arg[0] }, // data
        deseq2_args.map { arg -> arg[1] }, // formula
        deseq2_args.map { arg -> arg[2] }  // ref
    )

    // R: 10-dream ------------------------------------------------------------
    dream_formula = channel.fromList(params.dream_formula)
    dream_ref     = channel.fromList(params.dream_ref)
    dream_pairs   = dream_formula.merge(dream_ref)

    dream_args = pseudobulk.out.combine(dream_pairs)

    dream(
        dream_args.map { arg -> arg[0] }, // data
        dream_args.map { arg -> arg[1] }, // formula
        dream_args.map { arg -> arg[2] }  // ref
    )

    // DGE single cell ========================================================
    // R: 11-nebula -----------------------------------------------------------
    nebula_formula = channel.fromList(params.nebula_formula)
    nebula_ref     = channel.fromList(params.nebula_ref)
    nebula_pairs   = nebula_formula.merge(nebula_ref)

    nebula_args = subset.out.combine(nebula_pairs)

    nebula(
        nebula_args.map { arg -> arg[0] }, // data
        nebula_args.map { arg -> arg[1] }, // formula
        nebula_args.map { arg -> arg[2] }  // ref
    )

    // R: 12-sceptre ----------------------------------------------------------
    sceptre_formula = channel.fromList(params.sceptre_formula)
    sceptre_ref     = channel.fromList(params.sceptre_ref)
    sceptre_pairs   = sceptre_formula.merge(sceptre_ref)

    sceptre_args = subset.out
    .combine(sceptre_pairs)

    sceptre(
        sceptre_args.map { arg -> arg[0] }, // data
        sceptre_args.map { arg -> arg[1] }, // formula
        sceptre_args.map { arg -> arg[2] }  // ref
    )

    // R: 13-wilcox -----------------------------------------------------------
    wilcox_formula = channel.fromList(params.wilcox_formula)
    wilcox_ref     = channel.fromList(params.wilcox_ref)
    wilcox_pairs   = wilcox_formula.merge(wilcox_ref)

    wilcox_args = subset.out.combine(wilcox_pairs)

    wilcox(
        wilcox_args.map { arg -> arg[0] }, // data
        wilcox_args.map { arg -> arg[1] }, // formula
        wilcox_args.map { arg -> arg[2] }  // ref
    )

    // R: 14-scmageck ---------------------------------------------------------
    scmageck_formula = channel.fromList(params.scmageck_formula)
    scmageck_ref     = channel.fromList(params.scmageck_ref)
    scmageck_pairs   = scmageck_formula.merge(scmageck_ref)

    scmageck_args = subset.out.combine(scmageck_pairs)

    scmageck(
        scmageck_args.map { arg -> arg[0] }, // data
        scmageck_args.map { arg -> arg[1] }, // formula
        scmageck_args.map { arg -> arg[2] }  // ref
    )

    // Visualize ==============================================================
    // R: 15-visualize --------------------------------------------------------
    dge_results = limma.out.collect()
        .mix(edger.out)
        .mix(deseq2.out)
        .mix(dream.out)
        .mix(nebula.out)
        .mix(sceptre.out)
        .mix(wilcox.out)
        .mix(scmageck.out)
        .collect()

    visualize(
        dge_results
    )

    // Publish ================================================================
    publish:
    // Sim. res
    res_prepare                    = prepare.out
    res_splatter                   = splatter.out
    res_scdesign3_offtarget_unique = scdesign3_offtarget_unique.out
    res_compare                    = compare.out
    res_subset                     = subset.out
    res_pseudobulk                 = pseudobulk.out

    // DGE res
    res_limma    = limma.out
    res_edger    = edger.out
    res_deseq2   = deseq2.out
    res_dream    = dream.out
    res_nebula   = nebula.out
    res_sceptre  = sceptre.out
    res_wilcox   = wilcox.out
    res_scmageck = scmageck.out

    // Visualize
    res_visualize = visualize.out.panels
}

output {

    // Sim. res
    res_prepare {
        path "01-prepare"
        mode "copy"
    }
    res_splatter {
        path "02-splatter"
        mode "copy"
    }
    res_scdesign3_offtarget_unique {
        path "03-scdesign3_offtarget_unique"
        mode "copy"
    }
    res_compare {
        path "04-compare"
        mode "copy"
    }
    res_subset {
        path "05-subset"
        mode "copy"
    }
    res_pseudobulk {
        path "06-pseudobulk"
        mode "copy"
    }

    // DGE res
    res_limma {
        path "07-limma"
        mode "copy"
    }
    res_edger {
        path "08-edger"
        mode "copy"
    }
    res_deseq2 {
        path "09-deseq2"
        mode "copy"
    }
    res_dream {
        path "10-dream"
        mode "copy"
    }
    res_nebula {
        path "11-nebula"
        mode "copy"
    }
    res_sceptre {
        path "12-sceptre"
        mode "copy"
    }
    res_wilcox {
        path "13-wilcox"
        mode "copy"
    }
    res_scmageck {
        path "14-scmageck"
        mode "copy"
    }

    // Visualize
    res_visualize {
        path "15-visualize"
        mode "copy"
    }
}