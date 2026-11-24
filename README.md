

Differential gene expression (DGE) scripts:
The project benchmarked common DGE tools to analyse the performance on Perturb-Seq datasets.

#### 07-limma.R
Best practice limma-voom differential gene expression pipeline.

**Arguments:**
- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design. (e.g. `~ genotype + sample`).
- `ref`: The reference level for the target covariate. (e.g. `"NTC"`).

**Description:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Model matrix: `stats::model.matrix`, using the specified `formula`.
2. Gene filtering: `edgeR::filterByExpr`, using the the model matrix.
3. Normalization: `edgeR::calcNormFactors`, using library sizes of remaining genes.
4. Precision weights: `limma::voom`, using the model matrix.
5. Model fit: `limma::lmFit`.
6. Test: `limma::eBayes`.

NOTE: The differential statistics are merged with simulation ground truth data.

### **edgeR**

#### **Arguments**
* `formula`: An object of class `formula` specifying the design layout (e.g., `~ genotype`).
* `ref`: A string indicating the reference level for the target covariate (e.g., `"NTC"`).

#### **Description**
The **first covariate** specified in the `formula` is treated as the **target covariate**. The pipeline exclusively extracts differentially expressed genes (DEGs) associated with this target. 

To establish the baseline for comparison, the corresponding column in the metadata dataframe is converted into a factor. 
The level matching the `ref` argument is explicitly set as the **first level (reference baseline)**.

The analysis is executed through the following sequential workflow:
1. **Gene filtering**: Lowly expressed genes are filtered out using `edgeR::filterByExpr`, leveraging the structure of the model matrix.
2. **Normalization**:  Normalization factors are computed via `edgeR::calcNormFactors`, with library sizes recalculated based on the remaining genes.
3. **Dispersion**:     Empirical Bayes dispersion estimation using `edgeR::estimateDisp`.
4. **Model fit**:      Quasi-likelihood negative binomial generalized linear models are fitted to the counts using `edgeR::glmQLFit`.
5. **Test**:           Differential expression statistics are generated using `edgeR::glmQLFTest`.

Finally, the calculated differential statistics are merged with simulation ground truth data for downstream evaluation and publication.

### **DESeq2**

#### **Arguments**
* `formula`: An object of class `formula` specifying the design layout (e.g., `~ genotype`).
* `ref`: A string indicating the reference level for the target covariate (e.g., `"NTC"`).

#### **Description**
The **first covariate** specified in the `formula` is treated as the **target covariate**. The pipeline exclusively extracts differentially expressed genes (DEGs) associated with this target. 

To establish the baseline for comparison, the corresponding column in the metadata dataframe is converted into a factor. 
The level matching the `ref` argument is explicitly set as the **first level (reference baseline)**.

The analysis is executed through the following sequential workflow:
1. **Gene filtering**: Lowly expressed genes are filtered out using `edgeR::filterByExpr`, leveraging the structure of the model matrix.
2. **Dataset**:        A `DESeqDataSet` object is generated from the filtered raw count matrix and metadata via `DESeq2::DESeqDataSetFromMatrix`.
3. **Test**:           Differential expression statistics are generated using `DESeq2::DESeq`.

Finally, the calculated differential statistics are merged with simulation ground truth data for downstream evaluation and publication.

### **DREAM**

#### **Arguments**
* `formula`: An object of class `formula` specifying the design layout (e.g., `~ genotype`).
* `ref`: A string indicating the reference level for the target covariate (e.g., `"NTC"`).

#### **Description**
The **first covariate** specified in the `formula` is treated as the **target covariate**. The pipeline exclusively extracts differentially expressed genes (DEGs) associated with this target. 

To establish the baseline for comparison, the corresponding column in the metadata dataframe is converted into a factor. 
The level matching the `ref` argument is explicitly set as the **first level (reference baseline)**.

The analysis is executed through the following sequential workflow:
1. **Gene filtering**:    Lowly expressed genes are filtered out using `edgeR::filterByExpr`, leveraging the fixed-effects formulation of the model matrix.
3. **Normalization**:     Normalization factors are computed via `edgeR::calcNormFactors`, with library sizes recalculated based on the remaining genes.
3. **Precision weights**: Precision weights are calculated using `variancePartition::voomWithDreamWeights`.
4. **Model fit**:         Linear mixed models are fitted to the weights and expression data using `variancePartition::dream`.
5. **Test**:              Differential expression statistics are generated using `variancePartition::eBayes`.

Finally, the calculated differential statistics are merged with simulation ground truth data for downstream evaluation and publication.