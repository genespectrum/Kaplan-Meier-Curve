# ==========================================================
# Survival analysis pipeline using TCGA-READ RNA-seq data
# Gene: CDKN1A
# Samples analyzed: 130 primary tumors
# Author: <Your Name>
# ==========================================================

# --------------------------
# Libraries
# --------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(tibble)
library(tidyr)
library(survival)
library(survminer)

# --------------------------
# User-defined parameters
# --------------------------
TARGET_GENE   <- "CDKN1A"
N_SAMPLES     <- 130
TCGA_PROJECT  <- "TCGA-READ"

# --------------------------
# Function 1: Clinical data
# --------------------------
get_clinical_data <- function(project_id) {
  GDCquery_clinic(project_id, type = "clinical") %>%
    transmute(
      patient_id = submitter_id,
      status = vital_status,
      event = ifelse(vital_status == "Dead", 1, 0),
      time = ifelse(
        vital_status == "Dead",
        days_to_death,
        days_to_last_follow_up
      )
    )
}

clinical_df <- get_clinical_data(TCGA_PROJECT)

# --------------------------
# Function 2: Select samples
# --------------------------
select_primary_samples <- function(project_id, n) {
  q <- GDCquery(
    project = project_id,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor",
    access = "open"
  )
  
  getResults(q) %>%
    filter(sample_type == "Primary Tumor") %>%
    pull(cases) %>%
    unique() %>%
    sort() %>%
    head(n)
}

selected_samples <- select_primary_samples(TCGA_PROJECT, N_SAMPLES)

# --------------------------
# Function 3: Download & prep RNA-seq
# --------------------------
prepare_expression_data <- function(project_id, barcodes) {
  q <- GDCquery(
    project = project_id,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor",
    access = "open",
    barcode = barcodes
  )
  
  GDCdownload(q)
  GDCprepare(q)
}

se_object <- prepare_expression_data(TCGA_PROJECT, selected_samples)

# --------------------------
# Normalization
# --------------------------
counts_mat <- assay(se_object, "unstranded")
col_info   <- as.data.frame(colData(se_object))
gene_info  <- as.data.frame(rowData(se_object))

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = col_info,
  design    = ~ 1
)

dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = TRUE)
expr_vst <- assay(vsd)

# --------------------------
# Function 4: Gene extraction
# --------------------------
extract_gene_expression <- function(expr_matrix, gene_df, gene_symbol) {
  expr_matrix %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(
      cols = -gene_id,
      names_to = "sample_barcode",
      values_to = "expression"
    ) %>%
    left_join(gene_df, by = "gene_id") %>%
    filter(gene_name == gene_symbol)
}

gene_expr <- extract_gene_expression(expr_vst, gene_info, TARGET_GENE)

# --------------------------
# Grouping strategy
# --------------------------
gene_expr <- gene_expr %>%
  mutate(
    patient_id = sub("-01.*", "", sample_barcode),
    expr_group = ifelse(
      expression >= median(expression),
      "High expression",
      "Low expression"
    )
  )

# --------------------------
# Merge survival + expression
# --------------------------
survival_df <- gene_expr %>%
  inner_join(clinical_df, by = "patient_id")

# --------------------------
# Kaplan–Meier model
# --------------------------
km_fit <- survfit(
  Surv(time, event) ~ expr_group,
  data = survival_df
)

# --------------------------
# Visualization
# --------------------------
png(filename = "survival_plot_read_cdkn1a.png",
    width = 1200, height = 650)
ggsurvplot(
  km_fit,
  data = survival_df,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "CDKN1A",
  legend.labs = c("High expression", "Low expression"),
  title = "Kaplan–Meier Overall Survival Analysis of CDKN1A in TCGA-READ"
)
dev.off()
# --------------------------
# Log-rank statistics
# --------------------------
survdiff(
  Surv(time, event) ~ expr_group,
  data = survival_df
)


