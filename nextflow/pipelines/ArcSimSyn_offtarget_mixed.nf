// Modules
include { splatter }   from '../modules/splatter.nf'
include { scdesign3 }  from '../modules/scdesign3_offtarget_mixed.nf'
include { subset }     from '../modules/subset.nf'
include { pseudobulk } from '../modules/pseudobulk.nf'
include { limma }      from '../modules/limma.nf'
include { edger }      from '../modules/edger.nf'
include { deseq2 }     from '../modules/deseq2.nf'
include { dream }      from '../modules/dream.nf'
include { nebula }     from '../modules/nebula.nf'
include { sceptre }    from '../modules/sceptre.nf'
include { wilcox }     from '../modules/wilcox.nf'
include { scmageck }   from '../modules/scmageck.nf'

// Pipeline parameters
params {
    data: Path = "data/JB22_PM/rna_integrated.rds"
    
    // splatter
    splatter_facLoc = [ 0.1, 0.5, 1.5, 1, 2 ]
    
    // scdesign3
    scdesign3_nOn   = [ 10, 50, 100, 200, 500 ]
    scdesign3_koEff = [ 25, 50, 75, 100 ]

    // subset
    subset_nGuide = [ 10, 20, 50, 100, 200, 500, 1000 ]
    
    // limma
    limma_formula = [ '~genotype' , '~genotype + sample' , '~guide' , '~guide + sample' ]
    limma_ref     = [ 'NTC' , 'NTC' ,  'NTC-1' , 'NTC-1' ]
    
    // edger
    edger_formula = [ '~genotype' , '~genotype + sample' , '~guide' , '~guide + sample' ]
    edger_ref     = [ 'NTC' , 'NTC' ,  'NTC-1' , 'NTC-1' ]

    // deseq2
    deseq2_formula = [ '~genotype' , '~genotype + sample' , '~guide' , '~guide + sample' ]
    deseq2_ref     = [ 'NTC' , 'NTC' ,  'NTC-1' , 'NTC-1' ]
    
    // dream
    dream_formula = [ '~genotype' , '~genotype + sample' , '~genotype + (1|sample)' , '~genotype + sample + (1|guide)', '~genotype + (1|sample) + (1|guide)' , '~guide' , '~guide + sample' , '~guide + (1|sample)' ]
    dream_ref     = [ 'NTC' , 'NTC' , 'NTC' , 'NTC', 'NTC' , 'NTC-1' , 'NTC-1' , 'NTC-1' ]
    
    // nebula
    nebula_formula = [ '~genotype' , '~genotype + sample' , '~genotype + (1|sample)' , '~genotype + sample + (1|guide)', '~guide' , '~guide + sample' , '~guide + (1|sample)' ]
    nebula_ref     = [ 'NTC' , 'NTC' , 'NTC' , 'NTC', 'NTC-1' , 'NTC-1' , 'NTC-1' ]
    
    // sceptre
    sceptre_formula = [ '~genotype' , '~genotype + sample' ]
    sceptre_ref     = [ 'NTC' , 'NTC' ]
    
    // wilcox
    wilcox_formula = [ '~genotype' ]
    wilcox_ref     = [ 'NTC' ]

    // scmageck
    scmageck_formula = [ '~genotype' ]
    scmageck_ref     = [ 'NTC' ]
}

workflow {
    main:

    // Simualtions ============================================================

    // R: splatter ------------------------------------------------------------
    splatter_facLoc_ch = channel.fromList(params.splatter_facLoc)

    splatter_ch = splatter_facLoc_ch
        .combine(channel.value(params.data))
        .map { facLoc, file -> [file, facLoc] }
        // .view { filePath, facLoc -> 
        //     String.format(
        //         "SPLATTER   [data: %-50s facLoc: %-5s]", 
        //         file(filePath).name, facLoc
        //     )
        // }

    splatter(
        splatter_ch.map { table -> table[0] }, // data
        splatter_ch.map { table -> table[1] }  // facLoc
    )

    // R: scDesign3 -----------------------------------------------------------
    scdesign3_nOn_ch   = channel.fromList(params.scdesign3_nOn)
    scdesign3_koEff_ch = channel.fromList(params.scdesign3_koEff)

    //scdesign3_pairs_ch = scdesign3_nOn_ch.merge(scdesign3_koEff_ch)

    scdesign3_ch = splatter.out
        .map { out -> [ out[1], out[0] ] } // [data, facLoc]
        .combine(scdesign3_nOn_ch)
        .combine(scdesign3_koEff_ch)
        // .view { filePath, facLoc, nOn, koEff -> 
        //     String.format(
        //         "SCDESIGN3  [data: %-50s facLoc: %-5s nOn: %-5s koEff: %-5s]", 
        //         file(filePath).name, facLoc, nOn, koEff
        //     )
        // }

    scdesign3(
        scdesign3_ch.map { table -> table[0] }, // data 
        scdesign3_ch.map { table -> table[1] }, // facLoc 
        scdesign3_ch.map { table -> table[2] }, // nOn
        scdesign3_ch.map { table -> table[3] }  // koEff
    )

    // R: subset --------------------------------------------------------------
    subset_nGuide_ch = channel.fromList(params.subset_nGuide)

    subset_ch = scdesign3.out
        .map { out -> [ out[3], out[0], out[1], out[2] ] } // [data, facLoc, nOn, koEff]
        .combine(subset_nGuide_ch)
        // .view { filePath, facLoc, nOn, koEff, nGuide -> 
        //     String.format(
        //         "SUBSET     [data: %-50s facLoc: %-5s nOn: %-5s koEff: %-5s nGuide: %-5s]", 
        //         file(filePath).name, facLoc, nOn, koEff, nGuide
        //     )
        // }

    subset(
        subset_ch.map { table -> table[0] }, // data 
        subset_ch.map { table -> table[1] }, // facLoc 
        subset_ch.map { table -> table[2] }, // nOn
        subset_ch.map { table -> table[3] }, // koEff
        subset_ch.map { table -> table[4] }  // nGuide
    )

    // R: pseudobulk ----------------------------------------------------------
    pseudobulk_ch = subset.out
        .map { out -> [ out[4], out[0], out[1], out[2], out[3] ] } // [data, facLoc, nOn, koEff, nGuide]
        // .view { filePath, facLoc, nOn, koEff, nGuide -> 
        //     String.format(
        //         "PSEUDOBULK [data: %-50s facLoc: %-5s nOn: %-5s koEff: %-5s nGuide: %-5s]", 
        //         file(filePath).name, facLoc, nOn, koEff, nGuide
        //     )
        // }

    pseudobulk(
        pseudobulk_ch.map { table -> table[0] }, // data
        pseudobulk_ch.map { table -> table[1] }, // facLoc
        pseudobulk_ch.map { table -> table[2] }, // nOn
        pseudobulk_ch.map { table -> table[3] }, // koEff
        pseudobulk_ch.map { table -> table[4] }  // nGuide
    )

    // DGE pseudobulk =========================================================

    // limma ------------------------------------------------------------------
    limma_formula_ch = channel.fromList(params.limma_formula)
    limma_ref_ch     = channel.fromList(params.limma_ref)

    limma_pairs_ch   = limma_formula_ch.merge(limma_ref_ch)

    dge_limma_ch = pseudobulk.out
        .map { out -> [ out[4], out[0], out[1], out[2], out[3] ] } // [data, facLoc, nOn, koEff, nGuide]
        .combine(limma_pairs_ch)
        // .view { filePath, facLoc, nOn, koEff, nGuide, formula, ref -> 
        //     String.format(
        //         "LIMMA      [data: %-50s facLoc: %-5s nOn: %-5s koEff: %-5s nGuide: %-5s formula: %-35s ref: %-8s ]", 
        //         file(filePath).name, facLoc, nOn, koEff, nGuide, formula, ref
        //     )
        // }

    limma(
        dge_limma_ch.map { table -> table[0] }, // data
        dge_limma_ch.map { table -> table[1] }, // facLoc
        dge_limma_ch.map { table -> table[2] }, // nOn
        dge_limma_ch.map { table -> table[3] }, // koEff
        dge_limma_ch.map { table -> table[4] }, // nGuide
        dge_limma_ch.map { table -> table[5] }, // formula
        dge_limma_ch.map { table -> table[6] }  // ref
    )

    // edger ------------------------------------------------------------------
    edger_formula_ch = channel.fromList(params.edger_formula)
    edger_ref_ch     = channel.fromList(params.edger_ref)

    edger_pairs_ch   = edger_formula_ch.merge(edger_ref_ch)

    dge_edger_ch = pseudobulk.out
        .map { out -> [ out[4], out[0], out[1], out[2], out[3] ] } // [data, facLoc, nOn, koEff, nGuide]
        .combine(edger_pairs_ch)
        // .view { filePath, facLoc, nOn, koEff, nGuide, formula, ref -> 
        //     String.format(
        //         "EDGER      [data: %-50s facLoc: %-5s nOn: %-5s koEff: %-5s nGuide: %-5s formula: %-35s ref: %-8s ]", 
        //         file(filePath).name, facLoc, nOn, koEff, nGuide, formula, ref
        //     )
        // }

    edger(
        dge_edger_ch.map { table -> table[0] }, // data
        dge_edger_ch.map { table -> table[1] }, // facLoc
        dge_edger_ch.map { table -> table[2] }, // nOn
        dge_edger_ch.map { table -> table[3] }, // koEff
        dge_edger_ch.map { table -> table[4] }, // nGuide
        dge_edger_ch.map { table -> table[5] }, // formula
        dge_edger_ch.map { table -> table[6] }  // ref
    )


    // deseq2 -----------------------------------------------------------------
    deseq2_formula_ch = channel.fromList(params.deseq2_formula)
    deseq2_ref_ch     = channel.fromList(params.deseq2_ref)

    deseq2_pairs_ch   = deseq2_formula_ch.merge(deseq2_ref_ch)

    dge_deseq2_ch = pseudobulk.out
        .map { out -> [ out[4], out[0], out[1], out[2], out[3] ] } // [data, facLoc, nOn, koEff, nGuide]
        .combine(deseq2_pairs_ch)
        // .view { filePath, facLoc, nOn, koEff, nGuide, formula, ref -> 
        //     String.format(
        //         "DESEQ2     [data: %-50s facLoc: %-5s nOn: %-5s koEff: %-5s nGuide: %-5s formula: %-35s ref: %-8s ]", 
        //         file(filePath).name, facLoc, nOn, koEff, nGuide, formula, ref
        //     )
        // }

    deseq2(
        dge_deseq2_ch.map { table -> table[0] }, // data
        dge_deseq2_ch.map { table -> table[1] }, // facLoc
        dge_deseq2_ch.map { table -> table[2] }, // nOn
        dge_deseq2_ch.map { table -> table[3] }, // koEff
        dge_deseq2_ch.map { table -> table[4] }, // nGuide 
        dge_deseq2_ch.map { table -> table[5] }, // formula
        dge_deseq2_ch.map { table -> table[6] }  // ref
    )

    // dream ------------------------------------------------------------------
    dream_formula_ch = channel.fromList(params.dream_formula)
    dream_ref_ch     = channel.fromList(params.dream_ref)

    dream_pairs_ch   = dream_formula_ch.merge(dream_ref_ch)

    dge_dream_ch = pseudobulk.out
        .map { out -> [ out[4], out[0], out[1], out[2], out[3] ] } // [data, facLoc, nOn, koEff, nGuide]
        .combine(dream_pairs_ch)
        // .view { filePath, facLoc, nOn, koEff, nGuide, formula, ref -> 
        //     String.format(
        //         "DREAM      [data: %-50s facLoc: %-5s nOn: %-5s koEff: %-5s nGuide: %-5s formula: %-35s ref: %-8s ]", 
        //         file(filePath).name, facLoc, nOn, koEff, nGuide, formula, ref
        //     )
        // }

    dream(
        dge_dream_ch.map { table -> table[0] }, // data
        dge_dream_ch.map { table -> table[1] }, // facLoc
        dge_dream_ch.map { table -> table[2] }, // nOn
        dge_dream_ch.map { table -> table[3] }, // koEff
        dge_dream_ch.map { table -> table[4] }, // nGuide 
        dge_dream_ch.map { table -> table[5] }, // formula
        dge_dream_ch.map { table -> table[6] }  // ref
    )

    // DGE single cell ========================================================

    // nebula -----------------------------------------------------------------
    nebula_formula_ch = channel.fromList(params.nebula_formula)
    nebula_ref_ch     = channel.fromList(params.nebula_ref)

    nebula_pairs_ch   = nebula_formula_ch.merge(nebula_ref_ch)

    dge_nebula_ch = subset.out
        .map { out -> [ out[4], out[0], out[1], out[2], out[3] ] } // [data, facLoc, nOn, koEff, nGuide]
        .combine(nebula_pairs_ch)
        // .view { filePath, facLoc, nOn, koEff, nGuide, formula, ref -> 
        //     String.format(
        //         "NEBULA     [data: %-50s facLoc: %-5s nOn: %-5s koEff: %-5s nGuide: %-5s formula: %-35s ref: %-8s ]", 
        //         file(filePath).name, facLoc, nOn, koEff, nGuide, formula, ref
        //     )
        // }

    nebula(
        dge_nebula_ch.map { table -> table[0] }, // data
        dge_nebula_ch.map { table -> table[1] }, // facLoc
        dge_nebula_ch.map { table -> table[2] }, // nOn 
        dge_nebula_ch.map { table -> table[3] }, // koEff
        dge_nebula_ch.map { table -> table[4] }, // nGuide
        dge_nebula_ch.map { table -> table[5] }, // formula
        dge_nebula_ch.map { table -> table[6] }  // ref
    )

    // sceptre ----------------------------------------------------------------
    sceptre_formula_ch = channel.fromList(params.sceptre_formula)
    sceptre_ref_ch     = channel.fromList(params.sceptre_ref)

    sceptre_pairs_ch   = sceptre_formula_ch.merge(sceptre_ref_ch)

    dge_sceptre_ch = subset.out
        .map { out -> [ out[4], out[0], out[1], out[2], out[3] ] } // [data, facLoc, nOn, koEff, nGuide]
        .combine(sceptre_pairs_ch)
        // .view { filePath, facLoc, nOn, koEff, nGuide, formula, ref -> 
        //     String.format(
        //         "SCEPTRE    [data: %-50s facLoc: %-5s nOn: %-5s koEff: %-5s formula: %-35s ref: %-8s ]", 
        //         file(filePath).name, facLoc, nOn, koEff, nGuide, formula, ref
        //     )
        // }

    sceptre(
        dge_sceptre_ch.map { table -> table[0] }, // data
        dge_sceptre_ch.map { table -> table[1] }, // facLoc
        dge_sceptre_ch.map { table -> table[2] }, // nOn 
        dge_sceptre_ch.map { table -> table[3] }, // koEff
        dge_sceptre_ch.map { table -> table[4] }, // nGuide 
        dge_sceptre_ch.map { table -> table[5] }, // formula
        dge_sceptre_ch.map { table -> table[6] }  // ref
    )

    // wilcox -----------------------------------------------------------------
    wilcox_formula_ch = channel.fromList(params.wilcox_formula)
    wilcox_ref_ch     = channel.fromList(params.wilcox_ref)

    wilcox_pairs_ch   = wilcox_formula_ch.merge(wilcox_ref_ch)

    dge_wilcox_ch = subset.out
        .map { out -> [ out[4], out[0], out[1], out[2], out[3] ] } // [data, facLoc, nOn, koEff, nGuide]
        .combine(wilcox_pairs_ch)
        // .view { filePath, facLoc, nOn, koEff, nGuide, formula, ref -> 
        //     String.format(
        //         "WILCOX     [data: %-50s facLoc: %-5s nOn: %-5s koEff: %-5s nGuide: %-5s formula: %-35s ref: %-8s ]", 
        //         file(filePath).name, facLoc, nOn, koEff, nGuide, formula, ref
        //     )
        // }

    wilcox(
        dge_wilcox_ch.map { table -> table[0] }, // data
        dge_wilcox_ch.map { table -> table[1] }, // facLoc
        dge_wilcox_ch.map { table -> table[2] }, // nOn 
        dge_wilcox_ch.map { table -> table[3] }, // koEff
        dge_wilcox_ch.map { table -> table[4] }, // nGuide 
        dge_wilcox_ch.map { table -> table[5] }, // formula
        dge_wilcox_ch.map { table -> table[6] }  // ref
    )

    // scmageck ---------------------------------------------------------------
    scmageck_formula_ch = channel.fromList(params.scmageck_formula)
    scmageck_ref_ch     = channel.fromList(params.scmageck_ref)

    scmageck_pairs_ch   = scmageck_formula_ch.merge(scmageck_ref_ch)

    dge_scmageck_ch = subset.out
        .map { out -> [ out[4], out[0], out[1], out[2], out[3] ] } // [data, facLoc, nOn, koEff, nGuide]
        .combine(scmageck_pairs_ch)
        // .view { filePath, facLoc, nOn, koEff, nGuide, formula, ref -> 
        //     String.format(
        //         "SCMAGECK   [data: %-50s facLoc: %-5s nOn: %-5s koEff: %-5s nGuide: %-5s formula: %-35s ref: %-8s ]", 
        //         file(filePath).name, facLoc, nOn, koEff, nGuide, formula, ref
        //     )
        // }

    scmageck(
        dge_scmageck_ch.map { table -> table[0] }, // data
        dge_scmageck_ch.map { table -> table[1] }, // facLoc
        dge_scmageck_ch.map { table -> table[2] }, // nOn 
        dge_scmageck_ch.map { table -> table[3] }, // koEff
        dge_scmageck_ch.map { table -> table[4] }, // nGuide 
        dge_scmageck_ch.map { table -> table[5] }, // formula
        dge_scmageck_ch.map { table -> table[6] }  // ref
    )

    publish:
    // Sim. res
    res_splatter   = splatter.out
    res_scdesign3  = scdesign3.out
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
        path "nextflow/splatter"
        mode "copy"
    }
    res_scdesign3 {
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