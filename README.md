# Nextflow pipelines

Currently the repository holds 4 nextflow pipelines:
- ArcSimSyn_offtarget_none.nf
- ArcSimSyn_offtarget_unique.nf
- ArcSimSyn_offtarget_mixed.nf
- ArcSimSyn_singlecell_baseline.nf

## Run a nextflow pipeline
The pipeline execution is managed through bash scripts in the `bash` directory.
``` 
bash/ 
  ├── nxf_run_pipeline.sh    # Pipeline launch
  ├── nxf_submit_pipline.sh  # Slurm execution
  └── nxf_clean_pipeline.sh  # Cleans a pipeline
```
```
# Submit a new pipeline run
bash bash/nxf_run_pipeline.sh -pipeline ArcSimSyn_offtarget_none -run_id test001
```
Note: The current date is added to the `run_id` in the format `_YYYY-MM-DD`, e.g. `test001_2026-11-18`.
```
# Resume an existing run
bash bash/nxf_run_pipeline.sh -pipeline ArcSimSyn_offtarget_none -run_id test001_2026-11-18 -resume
```

**Output directories:**

Each pipeline execution is organized by pipeline name and run ID. This ensures that different workflows and individual executions remain separated. The following directory structure is used:
```
.nextflow/
└── PIPELINE_NAME/
    └── RUN_ID/
        └── cache/
            └── ...
            
work/
└── PIPELINE_NAME/
    └── RUN_ID/
        └── ...

results/
└── nextflow/
    └── PIPELINE_NAME/
        └── RUN_ID/
            └── ...

slurm/
└── PIPELINE_NAME/
    └── RUN_ID.log
```

## Pipeline execution workflow
The pipeline launcher `nxf_run_pipeline.sh` performs the following steps:
1. Pipeline selection: The user specifies a pipeline name and a custom run ID.
2. Run ID generation: For new runs, the current date (YYYY-MM-DD) is automatically appended.
3. SLURM submission: Submits the job using sbatch.

The submitted job is handled by `nxf_sub_pipeline.sh` via the following:

**I. New execution**
1. Checks that no existing Nextflow run directory exists for the selected pipeline and run ID.
2. Creates a dedicated Nextflow cache directory.
3. Creates a dedicated Netxflow work directory.
4. Creates a dedicated SLURM log file.

**II. Resume execution**
1. Checks that the original Nextflow cache directory exists.
2. Verifies that the cache contains at least one session directory.
3. Checks that the original working directory is available.
4. Reuses the existing cache and work directory to continue the execution.

Note: If any required directory or file is missing, the submission script terminates.

## Pipeline modules
### Differential expression methods

#### 07-limma.R

**Arguments:**

- `data`: The dataset in a SingleCellExperiment object.
- `formula`: The model formula specifying the design.
- `ref`: The reference level for the target covariate.

**Workflow:**

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

**Workflow:**

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

**Workflow:**

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

**Workflow:**

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


**Workflow:**

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


**Workflow:**

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


**Workflow:**

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


**Workflow:**

The first covariate specified in `formula` is treated as the *target covariate*. The pipeline exclusively extracts differentially expressed genes associated with this covariate. 

The corresponding column to the *target covariate* is converted into a factor. 
The value matching the `ref` argument is explicitly set as the *refernce*.

The analysis is executed through the following sequential workflow:
1. Guide assignment: A pseudo guide matrix is created.
1. Test: `scMAGeCK::scmageck_lr`.


NOTE: The differential statistics are merged with simulation ground truth data.