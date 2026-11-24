### Differential expression methods

#### 07-limma.R ---------------------------------------------------------------

**Arguments:**

- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design.
- `ref`: The reference level for the target covariate.

**Script:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Model matrix: `stats::model.matrix` using the specified formula.
2. Gene filtering: `edgeR::filterByExpr` using the the model matrix.
3. Normalization: `edgeR::calcNormFactors` using library sizes of remaining genes.
4. Precision weights: `limma::voom` using the model matrix.
5. Model fit: `limma::lmFit`.
6. Test: `limma::eBayes`.

NOTE: The differential statistics are merged with simulation ground truth data.

#### 08-edgeR.R ---------------------------------------------------------------

**Arguments:**

- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design.
- `ref`: The reference level for the target covariate.

**Script:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Model matrix: `stats::model.matrix` using the specified formula.
2. Gene filtering: `edgeR::filterByExpr` using the the model matrix.
3. Normalization: `edgeR::calcNormFactors` using library sizes of remaining genes.
4. Dispersion: `edgeR::estimateDisp` using the model matrix.
5. Model fit: `edgeR::glmQLFit` using the model matrix.
6. Test: `edgeR::glmQLFTest`.

#### 09-DESeq2.R --------------------------------------------------------------

**Arguments:**

- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design.
- `ref`: The reference level for the target covariate.

**Script:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Model matrix: `stats::model.matrix` using the specified formula.
2. Gene filtering: `edgeR::filterByExpr` using the the model matrix.
3. Dataset: `DESeq2::DESeqDataSetFromMatrix` using the specified formula.
4. Normalization: `DESeq2::estimateSizeFactors` using type = "poscounts".
5. Test: `DESeq2::DESeq`.

#### 10-dream-R ---------------------------------------------------------------

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
3. Precision weights: `variancePartition::voomWithDreamWeights` using the full formula.
4. **Model fit**:         Linear mixed models are fitted to the weights and expression data using `variancePartition::dream`.
5. **Test**:              Differential expression statistics are generated using `variancePartition::eBayes`.

Finally, the calculated differential statistics are merged with simulation ground truth data for downstream evaluation and publication.