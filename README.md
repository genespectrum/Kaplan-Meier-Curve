# Survival Analysis Plotting
### This repository provides an R-based pipeline that integrates TCGA-READ RNA-seq and clinical data to assess the prognostic significance of CDKN1A expression using Kaplan–Meier overall survival analysis.
## Overview
Survival analysis is a key approach in cancer genomics to understand whether gene expression levels are associated with patient outcomes. In this pipeline:
- TCGA-READ RNA-seq count data are downloaded and processed
- Gene expression is normalized using DESeq2 variance stabilizing transformation (VST)
- Patients are stratified into high and low expression groups based on median CDKN1A expression
- Overall survival is analyzed using Kaplan–Meier curves and log-rank tests

This script is designed to be reproducible, modular, and easy to adapt for other genes or TCGA projects.

### Dataset Details
- Project: TCGA-READ (Rectum Adenocarcinoma)
- Data type: RNA-seq (STAR – Counts)
- Sample type: Primary Tumor
- Number of samples analyzed: 130
- Target gene: CDKN1A

Clinical data are retrieved directly from the GDC and include:
- Vital status (Alive / Dead)
- Days to death
- Days to last follow-up

### Requirements
- R Packages
- Ensure the following R packages are installed:
- TCGAbiolinks
- SummarizedExperiment
- DESeq2
- dplyr
- tibble
- tidyr
- survival
- survminer

You can install missing packages using:

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2"))
install.packages(c("dplyr", "tibble", "tidyr", "survival", "survminer"))

### Step-by-Step Tutorial

1. Clinical Data Extraction

   Clinical metadata are downloaded from TCGA and processed to generate:

   event: survival event (1 = dead, 0 = alive)

   time: overall survival time in days

   This step ensures compatibility with survival analysis functions in R.

2. Sample Selection

   Only primary tumor samples are selected

   A fixed number of samples (130) is chosen for consistency

   Unique TCGA case barcodes are extracted

3. RNA-seq Download and Preparation

   RNA-seq raw count data are downloaded using GDCdownload()

   Data are converted into a SummarizedExperiment object using GDCprepare()

   This structure allows seamless integration of expression data and metadata.

4. Normalization

   Low-count genes are filtered out

   Variance Stabilizing Transformation (VST) is applied using DESeq2

   VST reduces technical noise and makes expression values comparable across samples.

5. Gene Expression Extraction

   Expression values for CDKN1A are extracted from the normalized matrix

   Data are reshaped into a long format

   Gene IDs are mapped to gene symbols

6. Expression-Based Grouping

   Patients are divided into two groups:

   High expression: expression ≥ median

   Low expression: expression < median

   This median-based stratification is commonly used in survival studies.

7. Survival Analysis

   Kaplan–Meier survival curves are generated

   Statistical significance is assessed using the log-rank test

   Risk tables display the number of patients at risk over time

### Output Files

survival_plot_read_cdkn1a.png

Kaplan–Meier overall survival plot showing high vs low CDKN1A expression groups, including:

- p-value
- risk table
- survival probability over time

### How to Run

1. Clone this repository

2. Open the R script in RStudio

3. Ensure internet access (required for TCGA download)

4. Run the script sequentially from top to bottom

5. Downloaded TCGA data will be stored automatically in your working directory.

### Customization

You can easily adapt this pipeline by modifying:

- TARGET_GENE – analyze a different gene
- N_SAMPLES – change the number of samples
- TCGA_PROJECT – use another TCGA cancer type (e.g., TCGA-BRCA)

### Use Cases
- Biomarker discovery
- Prognostic gene evaluation
- Survival analysis practice using TCGA data
