## Nextflow pipelines

Currently the repository holds 4 nextflow pipelines.
- ArcSimSyn_offtarget_none.nf
- ArcSimSyn_offtarget_unique.nf
- ArcSimSyn_offtarget_mixed.nf

- ArcSimSyn_singlecell_baseline.nf


### Differential expression methods

#### 07-limma.R

**Arguments:**

- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design.
- `ref`: The reference level for the target covariate.

**Script:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Model matrix: `model.matrix` using the specified formula.
2. Gene filtering: `edgeR::filterByExpr` using the the model matrix.
3. Normalization: `edgeR::calcNormFactors` using library sizes of remaining genes.
4. Precision weights: `limma::voom` using the model matrix.
5. Model fit: `limma::lmFit`.
6. Test: `limma::eBayes`.

#### 08-edgeR.R

**Arguments:**

- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design.
- `ref`: The reference level for the target covariate.

**Script:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Model matrix: `model.matrix` using the specified formula.
2. Gene filtering: `edgeR::filterByExpr` using the the model matrix.
3. Normalization: `edgeR::calcNormFactors` using library sizes of remaining genes.
4. Dispersion: `edgeR::estimateDisp` using the model matrix.
5. Model fit: `edgeR::glmQLFit` using the model matrix.
6. Test: `edgeR::glmQLFTest`.

#### 09-DESeq2.R

**Arguments:**

- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design.
- `ref`: The reference level for the target covariate.

**Script:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Model matrix: `model.matrix` using the specified formula.
2. Gene filtering: `edgeR::filterByExpr` using the the model matrix.
3. Dataset: `DESeq2::DESeqDataSetFromMatrix` using the specified formula.
4. Normalization: `DESeq2::estimateSizeFactors` using type = "poscounts".
5. Test: `DESeq2::DESeq`.

#### 10-dream-R

**Arguments:**

- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design.
- `ref`: The reference level for the target covariate.

**Script:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Formula: The specified formula is split into fixed effects and random effects.
2. Gene filtering: `edgeR::filterByExpr` using the fixed-effects model matrix.
3. Normalization: `edgeR::calcNormFactors` using library sizes of remaining genes.
4. Precision weights: `variancePartition::voomWithDreamWeights` using the full formula.
5. Model fit: `variancePartition::dream`.
6. Test: `variancePartition::eBayes`.

#### 11-nebula.R

**Arguments:**

- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design.
- `ref`: The reference level for the target covariate.


**Script:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Formula: The specified formula is split into fixed effects and random effects.
2. Model fit: `nebula::nebula` with cells being grouped by random effect.

#### 12-sceptre.R

**Arguments:**

- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design.
- `ref`: The reference level for the target covariate.


**Script:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Formula: The specified formula is used to create a covaritate dataframe.
2. Guide asignment: A pseudo guide matrix is created.
3. Test: Sceptre pipeline.


#### 13-wilcox.R

**Arguments:**

- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design.
- `ref`: The reference level for the target covariate.


**Script:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Normalization: `Seurat::NormalizeData`.
2. Test: `Seurat::FindMarkers`.


#### 14-scmageck.R

**Arguments:**

- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design.
- `ref`: The reference level for the target covariate.


**Script:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Guide assignment: A pseudo guide matrix is created.
1. Test: `scMAGeCK::scmageck_lr`.


NOTE: The differential statistics are merged with simulation ground truth data.