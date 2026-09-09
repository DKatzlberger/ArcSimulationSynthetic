

Differential gene expression (DGE) scripts:
The project benchmarked common DGE tools to analyse the performance on Perturb-Seq datasets.

### **Limma Pipeline**

#### **Arguments**
* `formula`: An object of class `formula` specifying the design layout (e.g., `~ genotype + sample`).
* `ref`: A string indicating the reference level for the target covariate (e.g., `"NTC"`).

#### **Description**
The **first covariate** specified in the `formula` is treated as the **target covariate**. The pipeline exclusively extracts differentially expressed genes (DEGs) associated with this target. 

To establish the baseline for comparison, the corresponding column in the metadata dataframe—retrieved from a `SingleCellExperiment` object using `SingleCellExperiment::colData`—is converted into a factor. The level matching the `ref` argument is explicitly set as the **first level (reference baseline)**.

The analysis is executed through the following sequential workflow:
1. **Design Matrix**: A model design matrix is constructed using the provided `formula`.
2. **Gene Filtering**: Lowly expressed genes are filtered out using `edgeR::filterByExpr`, leveraging the structure of the model matrix.
3. **Normalization**: Normalization factors are computed via `edgeR::calcNormFactors`, with library sizes automatically recalculated based on the remaining filtered genes.
4. **Precision Weights**: Observational precision weights are calculated using `limma::voom` to unlock RNA-seq data for linear modeling.
5. **Model Fitting**: Linear models are fitted to the weights and expression data using `limma::lmFit`.
6. **Empirical Bayes**: Differential expression statistics are smoothed and finalized by applying `limma::eBayes`.

Finally, the calculated differential statistics are merged with simulation ground truth data for downstream evaluation and publication.
