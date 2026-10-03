# =============================================================================
# ESCC Five-Gene Classifier: Independent Evaluation — Complete Analysis Script
#
# Study: Prognostic Evaluation and Multi-Cohort Cross-Platform Diagnostic
#        Replication of a Five-Gene Transcriptomic Classifier in Esophageal
#        Squamous Cell Carcinoma
#
# Authors: Muhammad Anas Akbar, Nawera Khan
# Affiliation: Department of Biotechnology, Government College University,
#              Lahore, Pakistan
# Contact: rnanasje@gmail.com
#
# Submitted to: PLOS ONE, September 2026
#
# Description:
#   This script reproduces all analyses reported in the manuscript.
#   It evaluates the published five-gene ESCC classifier
#   (SIM2, RFC4, COL1A1, MMP1, CST1; Khalil et al., Bioscience Insights
#   2026;3(2):478-495) across three publicly available cohorts:
#     1. TCGA-ESCA (RNA-seq; prognostic + exploratory diagnostic)
#     2. GSE20347 (Affymetrix GPL571; diagnostic replication, paired)
#     3. GSE38129 (Affymetrix GPL571; independent diagnostic replication, paired)
#
# R version: 4.6.1
# Required packages: GEOquery, survival, survminer, pROC
#
# Data availability:
#   TCGA-ESCA: https://xenabrowser.net (UCSC Xena platform)
#   GSE20347:  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20347
#   GSE38129:  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38129
#
# Run time: approximately 5-10 minutes (excluding TCGA download)
# =============================================================================

# =============================================================================
# 0. SETUP
# =============================================================================

# Install packages if needed (uncomment):
# install.packages(c("survival", "survminer", "pROC"))
# BiocManager::install("GEOquery")

library(GEOquery)
library(survival)
library(survminer)
library(pROC)

# Published coefficients from Khalil et al. (2026), Table 3
# Bioscience Insights, Vol 3(2), pp. 478-495
published_coefficients <- c(
  SIM2   = -0.824,
  RFC4   =  1.079,
  COL1A1 =  0.639,
  MMP1   =  0.712,
  CST1   =  0.521
)

genes_of_interest <- c("SIM2", "RFC4", "COL1A1", "MMP1", "CST1")

# GPL571 probe IDs (best probe per gene by highest mean expression)
final_probes_gpl571 <- c(
  SIM2   = "206558_at",
  RFC4   = "204023_at",
  COL1A1 = "202310_s_at",
  MMP1   = "204475_at",
  CST1   = "206224_at"
)

# =============================================================================
# 1. TCGA-ESCA DATA
# =============================================================================
# Data downloaded from UCSC Xena platform (https://xenabrowser.net)
# Dataset: TCGA Esophageal Carcinoma (TCGA-ESCA)
# Files required:
#   - Gene expression: TCGA-ESCA.htseq_fpkm-uq.tsv
#   - Clinical: Survival_SupplementalTable_S1_20171025_xena_sp.tsv
#   - Histology: TCGA-ESCA.GDC_phenotype.tsv
#
# NOTE: TCGA files must be downloaded manually from Xena before running
#       this section. GDC manifest (MANIFEST.txt) is provided in this
#       repository for the raw RNA-seq files.

# --- Load TCGA data (adjust paths to your downloaded files) ------------------
# tcga_expr     <- read.table("TCGA-ESCA.htseq_fpkm-uq.tsv",
#                              header=TRUE, sep="\t", row.names=1,
#                              check.names=FALSE)
# tcga_clinical <- read.table("Survival_SupplementalTable_S1.tsv",
#                              header=TRUE, sep="\t")
# tcga_pheno    <- read.table("TCGA-ESCA.GDC_phenotype.tsv",
#                              header=TRUE, sep="\t")

# --- Filter to ESCC histology ------------------------------------------------
# Step 1: identify ESCC samples from phenotype annotation
# escc_samples <- tcga_pheno$submitter_id.samples[
#   grepl("squamous", tcga_pheno$histological_type, ignore.case = TRUE)
# ]
# cat("ESCC samples in phenotype:", length(escc_samples), "\n")
# Expected: 100 squamous samples out of 204 total

# Step 2: extract five-gene expression for ESCC samples only
# escc_expr_5genes <- tcga_expr[
#   genes_of_interest,
#   intersect(escc_samples, colnames(tcga_expr)),
#   drop = FALSE
# ]

# Step 3: merge expression with clinical survival data
# merged_data_escc <- merge(
#   as.data.frame(t(escc_expr_5genes)),
#   tcga_clinical[, c("sample", "OS", "OS.time")],
#   by.x = "row.names", by.y = "sample"
# )
# colnames(merged_data_escc)[1] <- "sample"
# cat("After expression-survival merge:", nrow(merged_data_escc), "\n")
# Expected: 99 (one sample lost due to missing survival data)

# Step 4: diagnostic subset
# Barcode suffix: -01 = primary tumour, -11 = solid tissue normal
# -06 = metastatic (excluded)
# merged_data_diag <- as.data.frame(t(escc_expr_5genes))
# merged_data_diag$sample <- rownames(merged_data_diag)
# merged_data_diag$tissue_type <- ifelse(
#   grepl("-01$", merged_data_diag$sample), "Tumor",
#   ifelse(grepl("-11$", merged_data_diag$sample), "Normal", "Other")
# )
# merged_data_diag <- merged_data_diag[
#   merged_data_diag$tissue_type %in% c("Tumor", "Normal"), ]
# cat("TCGA diagnostic subset:\n")
# cat(table(merged_data_diag$tissue_type), "\n")
# Expected: Tumor=95, Normal=3

# Keep primary tumour only for survival analysis
surv_data <- merged_data_escc[
  grepl("-01$", merged_data_escc$sample),
]

# Ensure numeric expression
for (gene in genes_of_interest) {
  surv_data[[gene]] <- as.numeric(surv_data[[gene]])
}

# Within-cohort z-standardization (survival cohort)
surv_data$z_SIM2_surv   <- as.numeric(scale(surv_data$SIM2))
surv_data$z_RFC4_surv   <- as.numeric(scale(surv_data$RFC4))
surv_data$z_COL1A1_surv <- as.numeric(scale(surv_data$COL1A1))
surv_data$z_MMP1_surv   <- as.numeric(scale(surv_data$MMP1))
surv_data$z_CST1_surv   <- as.numeric(scale(surv_data$CST1))

# Apply published coefficients without refitting
surv_data$frozen_score_surv <-
  (published_coefficients["SIM2"]   * surv_data$z_SIM2_surv)   +
  (published_coefficients["RFC4"]   * surv_data$z_RFC4_surv)   +
  (published_coefficients["COL1A1"] * surv_data$z_COL1A1_surv) +
  (published_coefficients["MMP1"]   * surv_data$z_MMP1_surv)   +
  (published_coefficients["CST1"]   * surv_data$z_CST1_surv)

# Keep complete cases only
surv_data <- surv_data[
  complete.cases(
    surv_data$OS.time, surv_data$OS,
    surv_data$SIM2, surv_data$RFC4,
    surv_data$COL1A1, surv_data$MMP1,
    surv_data$CST1, surv_data$frozen_score_surv
  ),
]

cat("Survival cohort: n =", nrow(surv_data), "\n")
cat("OS events:", sum(surv_data$OS), "\n")
# Expected: n = 95, events = 32

# =============================================================================
# 2. SURVIVAL ANALYSIS — TCGA ESCC
# =============================================================================

# --- Univariable Cox: individual genes ---------------------------------------
cox_sim2   <- coxph(Surv(OS.time, OS) ~ SIM2,   data = surv_data)
cox_rfc4   <- coxph(Surv(OS.time, OS) ~ RFC4,   data = surv_data)
cox_col1a1 <- coxph(Surv(OS.time, OS) ~ COL1A1, data = surv_data)
cox_mmp1   <- coxph(Surv(OS.time, OS) ~ MMP1,   data = surv_data)
cox_cst1   <- coxph(Surv(OS.time, OS) ~ CST1,   data = surv_data)

cat("\n=== Table 2: Univariable Cox — Individual Genes ===\n")
# Expected: SIM2 HR=0.885 p=0.516 C=0.565
#           RFC4 HR=0.947 p=0.801 C=0.466
#           COL1A1 HR=0.931 p=0.481 C=0.570
#           MMP1 HR=0.932 p=0.399 C=0.540
#           CST1 HR=0.971 p=0.678 C=0.585
cox_gene_list <- list(
  SIM2=cox_sim2, RFC4=cox_rfc4, COL1A1=cox_col1a1,
  MMP1=cox_mmp1, CST1=cox_cst1
)
table2_results <- lapply(names(cox_gene_list), function(gene) {
  obj <- cox_gene_list[[gene]]
  s   <- summary(obj)
  data.frame(
    Gene        = gene,
    HR          = round(s$coefficients[,"exp(coef)"], 3),
    CI_lower    = round(s$conf.int[,"lower .95"], 3),
    CI_upper    = round(s$conf.int[,"upper .95"], 3),
    p_value     = round(s$coefficients[,"Pr(>|z|)"], 3),
    Concordance = round(s$concordance["C"], 3),
    stringsAsFactors = FALSE
  )
})
table2_results <- do.call(rbind, table2_results)
print(table2_results)
write.csv(table2_results, "Table2_univariable_Cox_results.csv", row.names=FALSE)

# --- Multivariable Cox: all five genes simultaneously ------------------------
cox_combined <- coxph(
  Surv(OS.time, OS) ~ SIM2 + RFC4 + COL1A1 + MMP1 + CST1,
  data = surv_data
)
cat("\n=== Table 4: Multivariable Cox (all 5 genes) ===\n")
# Expected: LR p=0.855, Concordance=0.586
s_combined <- summary(cox_combined)
cat("LR p-value:", round(s_combined$logtest["pvalue"], 3), "\n")
cat("Concordance:", round(s_combined$concordance["C"], 3), "\n")
print(s_combined$coefficients)

# Save Table 4 multivariable results (S2 Table in manuscript)
s2_table <- as.data.frame(s_combined$coefficients)
s2_table$HR        <- exp(s2_table$coef)
s2_table$CI_lower  <- s_combined$conf.int[,"lower .95"]
s2_table$CI_upper  <- s_combined$conf.int[,"upper .95"]
write.csv(s2_table, "S2_Table_multivariable_Cox.csv", row.names=TRUE)

# --- Coefficient-fixed score Cox ---------------------------------------------
cox_frozen <- coxph(
  Surv(OS.time, OS) ~ frozen_score_surv,
  data = surv_data
)
cat("\n=== Table 4: Coefficient-Fixed Score Cox ===\n")
# Expected: HR=0.962, 95% CI 0.817-1.132, p=0.641, C-index=0.546
s_frozen <- summary(cox_frozen)
cat("HR:", round(s_frozen$coefficients[,"exp(coef)"], 3), "\n")
cat("95% CI:", round(s_frozen$conf.int[,"lower .95"], 3), "-",
    round(s_frozen$conf.int[,"upper .95"], 3), "\n")
cat("p-value:", round(s_frozen$coefficients[,"Pr(>|z|)"], 3), "\n")
cat("C-index:", round(s_frozen$concordance["C"], 3), "\n")

# --- Proportional hazards test -----------------------------------------------
# Expected: frozen score p=0.53, combined model global p=0.38
ph_frozen   <- cox.zph(cox_frozen)
ph_combined <- cox.zph(cox_combined)
cat("\n=== PH Test (frozen score) global p =",
    round(ph_frozen$table["GLOBAL","p"], 3), "\n")
cat("=== PH Test (combined model) global p =",
    round(ph_combined$table["GLOBAL","p"], 3), "\n")
cat("Individual gene PH tests:\n")
print(round(ph_combined$table, 3))

# --- Pairwise Pearson correlations -------------------------------------------
# Expected: COL1A1-MMP1 r=0.58, SIM2-COL1A1 r=-0.41,
#           SIM2-MMP1 r=-0.37, COL1A1-CST1 r=0.42
cor_matrix <- cor(
  surv_data[, genes_of_interest],
  method = "pearson",
  use = "complete.obs"
)
cat("\n=== Pairwise Pearson Correlations ===\n")
print(round(cor_matrix, 2))
write.csv(as.data.frame(round(cor_matrix, 2)),
          "Pearson_correlation_matrix.csv", row.names=TRUE)

# --- Kaplan-Meier: frozen score (Fig 1) --------------------------------------
surv_data$frozen_group <- factor(
  ifelse(
    surv_data$frozen_score_surv >= median(surv_data$frozen_score_surv),
    "High frozen score", "Low frozen score"
  ),
  levels = c("Low frozen score", "High frozen score")
)

fit_frozen <- survfit(Surv(OS.time, OS) ~ frozen_group, data = surv_data)

# Log-rank test for frozen score KM (p=0.66 expected)
lr_frozen <- survdiff(Surv(OS.time, OS) ~ frozen_group, data = surv_data)
lr_frozen_p <- 1 - pchisq(lr_frozen$chisq, df = 1)
cat("\nFrozen score KM log-rank p =", round(lr_frozen_p, 3), "\n")
# Expected: p = 0.66

plot_frozen <- ggsurvplot(
  fit_frozen, data = surv_data,
  risk.table = TRUE, risk.table.height = 0.25,
  pval = TRUE, conf.int = FALSE,
  censor = TRUE, censor.shape = 3, censor.size = 3,
  xlab = "Overall survival time (days)",
  ylab = "Overall survival probability",
  legend.title = "Frozen five-gene score",
  legend.labs = c("Low", "High"),
  palette = c("#2166AC", "#B2182B"),
  ggtheme = theme_classic(),
  break.time.by = 500,
  title = "Overall survival by frozen five-gene classifier score"
)
png("Fig1_Frozen_Score_Kaplan_Meier.png",
    width = 2400, height = 2400, res = 300)
print(plot_frozen)
dev.off()

# --- Kaplan-Meier: individual genes (S1-S5 Fig) ------------------------------
make_gene_km <- function(data, gene_name, output_file) {
  group_col <- paste0(gene_name, "_group")
  data[[group_col]] <- factor(
    ifelse(
      data[[gene_name]] >= median(data[[gene_name]], na.rm = TRUE),
      paste("High", gene_name), paste("Low", gene_name)
    ),
    levels = c(paste("Low", gene_name), paste("High", gene_name))
  )
  fit <- survfit(
    as.formula(paste0("Surv(OS.time, OS) ~ ", group_col)),
    data = data
  )
  km_plot <- ggsurvplot(
    fit, data = data,
    risk.table = TRUE, risk.table.height = 0.25,
    pval = TRUE, conf.int = FALSE,
    censor = TRUE, censor.shape = 3, censor.size = 3,
    xlab = "Overall survival time (days)",
    ylab = "Overall survival probability",
    legend.title = paste(gene_name, "expression"),
    legend.labs = c("Low", "High"),
    palette = c("#2166AC", "#B2182B"),
    ggtheme = theme_classic(),
    break.time.by = 500,
    title = paste("Overall survival by", gene_name, "expression")
  )
  png(output_file, width = 2400, height = 2400, res = 300)
  print(km_plot)
  dev.off()
  return(km_plot)
}

plot_SIM2   <- make_gene_km(surv_data, "SIM2",
                             "S1_Fig_SIM2_Kaplan_Meier.png")
plot_RFC4   <- make_gene_km(surv_data, "RFC4",
                             "S2_Fig_RFC4_Kaplan_Meier.png")
plot_COL1A1 <- make_gene_km(surv_data, "COL1A1",
                             "S3_Fig_COL1A1_Kaplan_Meier.png")
plot_MMP1   <- make_gene_km(surv_data, "MMP1",
                             "S4_Fig_MMP1_Kaplan_Meier.png")
plot_CST1   <- make_gene_km(surv_data, "CST1",
                             "S5_Fig_CST1_Kaplan_Meier.png")

# --- Table 3: Log-rank p-values for individual gene KM plots -----------------
# Extract log-rank p-values from each gene KM analysis
table3_logrank <- lapply(genes_of_interest, function(gene) {
  group_col <- paste0(gene, "_group")
  surv_data[[group_col]] <- factor(
    ifelse(
      surv_data[[gene]] >= median(surv_data[[gene]], na.rm = TRUE),
      paste("High", gene), paste("Low", gene)
    ),
    levels = c(paste("Low", gene), paste("High", gene))
  )
  fit <- survfit(
    as.formula(paste0("Surv(OS.time, OS) ~ ", group_col)),
    data = surv_data
  )
  # Log-rank test
  lr_test <- survdiff(
    as.formula(paste0("Surv(OS.time, OS) ~ ", group_col)),
    data = surv_data
  )
  lr_p <- 1 - pchisq(lr_test$chisq, df = 1)

  # N high / low
  tab <- table(surv_data[[group_col]])
  events_tab <- tapply(surv_data$OS, surv_data[[group_col]], sum)

  data.frame(
    Gene = gene,
    N_High = unname(tab[paste("High", gene)]),
    N_Low  = unname(tab[paste("Low",  gene)]),
    Events_High = unname(events_tab[paste("High", gene)]),
    Events_Low  = unname(events_tab[paste("Low",  gene)]),
    Logrank_p = lr_p,
    stringsAsFactors = FALSE
  )
})
table3_logrank <- do.call(rbind, table3_logrank)
cat("\n=== Table 3: KM Log-rank p-values (individual genes) ===\n")
print(table3_logrank)
# Expected p-values:
# SIM2=0.3709, RFC4=0.7055, COL1A1=0.7590, MMP1=0.9108, CST1=0.4697

write.csv(table3_logrank,
          "Table3_KM_logrank_pvalues.csv", row.names = FALSE)

cat("Survival analysis complete.\n")

# =============================================================================
# 3. GEO DATA ACQUISITION
# =============================================================================

cat("\nDownloading GSE20347...\n")
gse20347      <- getGEO("GSE20347", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
gse20347_expr <- exprs(gse20347)
gse20347_pheno <- pData(gse20347)

# Tissue grouping from source_name_ch1
gse20347_pheno$tissue_group <- ifelse(
  grepl("tumor|tumour|carcinoma",
        gse20347_pheno$source_name_ch1, ignore.case = TRUE),
  "Tumor", "Normal"
)

# Patient IDs from sample title (e.g. E1507T -> E1507, E1507N -> E1507)
gse20347_pheno$patient_id <- sub(
  "[NT]$", "", gse20347_pheno$title
)

cat("GSE20347: n =", nrow(gse20347_pheno), "\n")
cat(table(gse20347_pheno$tissue_group), "\n")
# Expected: Normal=17, Tumor=17

cat("\nDownloading GSE38129...\n")
gse38129       <- getGEO("GSE38129", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
gse38129_expr  <- exprs(gse38129)
gse38129_pheno <- pData(gse38129)

# Tissue grouping from tissue type field
gse38129_pheno$tissue_group <- ifelse(
  grepl("tumor", gse38129_pheno$`tissue type:ch1`, ignore.case = TRUE),
  "Tumor", "Normal"
)

# Patient IDs from title (e.g. E10T -> E10, E10N -> E10)
gse38129_pheno$patient_id <- sub(
  "[NT]$", "", gse38129_pheno$title
)

cat("GSE38129: n =", nrow(gse38129_pheno), "\n")
cat(table(gse38129_pheno$tissue_group), "\n")
# Expected: Normal=30, Tumor=30

# =============================================================================
# 4. SCORE CALCULATION FUNCTION
# =============================================================================

compute_coef_fixed_score <- function(expr_matrix,
                                      probes,
                                      coefficients) {
  # Extract probe rows
  gene_expr <- expr_matrix[probes, , drop = FALSE]
  rownames(gene_expr) <- names(probes)
  # Within-cohort z-standardization per gene
  z_scores <- t(scale(t(gene_expr)))
  # Apply published coefficients (no refitting)
  score <- colSums(coefficients * z_scores)
  return(score)
}

# =============================================================================
# 5. GSE20347 DIAGNOSTIC ANALYSIS
# =============================================================================

# --- Compute coefficient-fixed score -----------------------------------------
gse20347_score <- compute_coef_fixed_score(
  gse20347_expr, final_probes_gpl571, published_coefficients
)

gse20347_merged <- data.frame(
  geo_accession = names(gse20347_score),
  coef_fixed_score = gse20347_score,
  tissue_group = gse20347_pheno[
    match(names(gse20347_score), rownames(gse20347_pheno)),
    "tissue_group"],
  patient_id = gse20347_pheno[
    match(names(gse20347_score), rownames(gse20347_pheno)),
    "patient_id"],
  stringsAsFactors = FALSE
)

# --- ROC analysis (Fig 3) ----------------------------------------------------
roc_gse20347 <- roc(
  response  = ifelse(gse20347_merged$tissue_group == "Tumor", 1, 0),
  predictor = gse20347_merged$coef_fixed_score,
  ci = TRUE, ci.method = "delong"
)
auc_gse20347 <- auc(roc_gse20347)
cat("\nGSE20347 AUC:", round(as.numeric(auc_gse20347), 4), "\n")
cat("DeLong 95% CI:", round(ci(roc_gse20347)[1], 4), "-",
    round(ci(roc_gse20347)[3], 4), "\n")
# Expected: AUC = 1.000

png("Fig3_GSE20347_ROC.png", width = 1800, height = 1500, res = 300)
plot(roc_gse20347, col = "#1B7837", lwd = 3,
     legacy.axes = TRUE, xlim = c(1, 0), ylim = c(0, 1),
     main = "GSE20347: Coefficient-Fixed, Cohort-Standardized Score")
abline(a = 0, b = 1, lty = 2, col = "grey70")
text(x = 0.55, y = 0.15, adj = 0, cex = 0.9,
     labels = paste0("AUC = 1.000\n17 ESCC / 17 adjacent normal\n",
                     "All 17 pairs: tumour score higher\n",
                     "Paired Wilcoxon p = 1.53 \u00d7 10\u207b\u2075"))
dev.off()

# --- Individual gene paired Wilcoxon tests (Table 7) ------------------------
gse20347_gene_expr <- gse20347_expr[final_probes_gpl571, , drop = FALSE]
rownames(gse20347_gene_expr) <- names(final_probes_gpl571)

gse20347_gene_long <- data.frame(
  sample_id = colnames(gse20347_gene_expr),
  t(gse20347_gene_expr),
  check.names = FALSE, stringsAsFactors = FALSE
)
gse20347_gene_long$tissue_group <- gse20347_pheno[
  match(gse20347_gene_long$sample_id, rownames(gse20347_pheno)),
  "tissue_group"]
gse20347_gene_long$patient_id <- gse20347_pheno[
  match(gse20347_gene_long$sample_id, rownames(gse20347_pheno)),
  "patient_id"]

paired_results_gse20347 <- lapply(genes_of_interest, function(gene) {
  normal_vals <- gse20347_gene_long[
    gse20347_gene_long$tissue_group == "Normal", c("patient_id", gene)]
  tumor_vals  <- gse20347_gene_long[
    gse20347_gene_long$tissue_group == "Tumor",  c("patient_id", gene)]
  colnames(normal_vals)[2] <- "Normal"
  colnames(tumor_vals)[2]  <- "Tumor"
  paired <- merge(normal_vals, tumor_vals, by = "patient_id")
  paired$diff <- paired$Tumor - paired$Normal
  test <- wilcox.test(paired$Tumor, paired$Normal,
                      paired = TRUE, exact = TRUE)
  data.frame(
    Gene = gene, N_pairs = nrow(paired),
    Normal_mean = mean(paired$Normal), Tumor_mean = mean(paired$Tumor),
    Mean_diff = mean(paired$diff), Median_diff = median(paired$diff),
    Pairs_concordant = sum(paired$diff > 0),
    Paired_Wilcoxon_p = test$p.value, stringsAsFactors = FALSE
  )
})
paired_results_gse20347 <- do.call(rbind, paired_results_gse20347)
paired_results_gse20347$BH_FDR <- p.adjust(
  paired_results_gse20347$Paired_Wilcoxon_p, method = "BH"
)
cat("\nGSE20347 Paired Gene Tests (Table 7):\n")
# Expected BH-FDR values:
# SIM2   raw p=7.63e-5, BH-FDR=9.54e-5  (15/17 pairs concordant)
# RFC4   raw p=1.53e-5, BH-FDR=3.81e-5  (17/17 pairs concordant)
# COL1A1 raw p=4.58e-5, BH-FDR=7.63e-5  (16/17 pairs concordant)
# MMP1   raw p=1.53e-5, BH-FDR=3.81e-5  (17/17 pairs concordant)
# CST1   raw p=1.07e-4, BH-FDR=1.07e-4  (15/17 pairs concordant)
print(paired_results_gse20347)

# --- Paired composite score test (GSE20347) ----------------------------------
paired_wide_gse20347 <- reshape(
  gse20347_merged[, c("patient_id", "tissue_group", "coef_fixed_score")],
  idvar = "patient_id", timevar = "tissue_group", direction = "wide"
)
paired_wide_gse20347$Tumor_minus_Normal <-
  paired_wide_gse20347$coef_fixed_score.Tumor -
  paired_wide_gse20347$coef_fixed_score.Normal

wilcox_score_gse20347 <- wilcox.test(
  paired_wide_gse20347$coef_fixed_score.Tumor,
  paired_wide_gse20347$coef_fixed_score.Normal,
  paired = TRUE, alternative = "two.sided", exact = TRUE
)
cat("\nGSE20347 paired score Wilcoxon p =",
    format(wilcox_score_gse20347$p.value, scientific = TRUE), "\n")
cat("Pairs with higher tumour score:",
    sum(paired_wide_gse20347$Tumor_minus_Normal > 0), "/ 17\n")
# Expected: p = 1.53e-5, 17/17

write.csv(paired_results_gse20347,
          "GSE20347_five_gene_paired_Wilcoxon_results.csv",
          row.names = FALSE)
write.csv(gse20347_merged,
          "GSE20347_coefficient_fixed_transfer_data.csv",
          row.names = FALSE)
write.csv(paired_wide_gse20347[, c("patient_id",
                                    "coef_fixed_score.Normal",
                                    "coef_fixed_score.Tumor",
                                    "Tumor_minus_Normal")],
          "GSE20347_signature_score_paired_differences.csv",
          row.names = FALSE)

# =============================================================================
# 6. GSE38129 DIAGNOSTIC ANALYSIS
# =============================================================================

# --- Compute coefficient-fixed score -----------------------------------------
gse38129_score <- compute_coef_fixed_score(
  gse38129_expr, final_probes_gpl571, published_coefficients
)

gse38129_merged <- data.frame(
  geo_accession = names(gse38129_score),
  title = gse38129_pheno[
    match(names(gse38129_score), rownames(gse38129_pheno)), "title"],
  coef_fixed_score = gse38129_score,
  tissue_group = gse38129_pheno[
    match(names(gse38129_score), rownames(gse38129_pheno)),
    "tissue_group"],
  patient_id = gse38129_pheno[
    match(names(gse38129_score), rownames(gse38129_pheno)),
    "patient_id"],
  stringsAsFactors = FALSE
)

# --- ROC analysis (Fig 4) ----------------------------------------------------
roc_gse38129 <- roc(
  response  = ifelse(gse38129_merged$tissue_group == "Tumor", 1, 0),
  predictor = gse38129_merged$coef_fixed_score,
  ci = TRUE, ci.method = "delong"
)
auc_gse38129    <- auc(roc_gse38129)
ci_delong_gse38129 <- ci(roc_gse38129)

set.seed(42)
ci_boot_gse38129 <- ci.auc(roc_gse38129, method = "bootstrap",
                             boot.n = 2000, stratified = TRUE)

cat("\nGSE38129 AUC:", round(as.numeric(auc_gse38129), 4), "\n")
cat("DeLong 95% CI:", round(ci_delong_gse38129[1], 4), "-",
    round(ci_delong_gse38129[3], 4), "\n")
cat("Bootstrap 95% CI:", round(ci_boot_gse38129[1], 4), "-",
    round(ci_boot_gse38129[3], 4), "\n")
# Expected: AUC=0.9611
# DeLong 95% CI: 0.9107-1.0000
# Bootstrap 95% CI: 0.9011-0.9978

png("Fig4_GSE38129_ROC.png", width = 1800, height = 1500, res = 300)
plot(roc_gse38129, col = "#1B7837", lwd = 3,
     legacy.axes = TRUE, xlim = c(1, 0), ylim = c(0, 1),
     main = "GSE38129: Coefficient-Fixed, Cohort-Standardized Score")
abline(a = 0, b = 1, lty = 2, col = "grey70")
text(x = 0.62, y = 0.18, adj = 0, cex = 0.95,
     labels = paste0(
       "AUC = ", format(round(as.numeric(auc_gse38129), 3), nsmall = 3),
       "\n30 ESCC / 30 adjacent normal"
     ))
dev.off()

# --- Unpaired Welch tests across 60 samples (Table 8) -----------------------
gse38129_gene_expr <- gse38129_expr[final_probes_gpl571, , drop = FALSE]
rownames(gse38129_gene_expr) <- names(final_probes_gpl571)

gse38129_gene_long <- data.frame(
  sample_id = colnames(gse38129_gene_expr),
  t(gse38129_gene_expr),
  check.names = FALSE, stringsAsFactors = FALSE
)
gse38129_gene_long$tissue_type <- gse38129_pheno[
  match(gse38129_gene_long$sample_id, rownames(gse38129_pheno)),
  "tissue_group"]
gse38129_gene_long$patient_id <- gse38129_pheno[
  match(gse38129_gene_long$sample_id, rownames(gse38129_pheno)),
  "patient_id"]

gene_tests_gse38129 <- lapply(genes_of_interest, function(gene) {
  test <- t.test(
    gse38129_gene_long[[gene]] ~ gse38129_gene_long$tissue_type
  )
  normal_mean <- unname(test$estimate["mean in group Normal"])
  tumor_mean  <- unname(test$estimate["mean in group Tumor"])
  data.frame(
    Gene = gene, Normal_mean = normal_mean, Tumor_mean = tumor_mean,
    Difference_Normal_minus_Tumor = normal_mean - tumor_mean,
    Direction_in_tumor = ifelse(tumor_mean > normal_mean, "Higher", "Lower"),
    Welch_p_value = test$p.value, stringsAsFactors = FALSE
  )
})
gene_tests_gse38129 <- do.call(rbind, gene_tests_gse38129)
gene_tests_gse38129$BH_FDR <- p.adjust(
  gene_tests_gse38129$Welch_p_value, method = "BH"
)
cat("\nGSE38129 Welch tests (Table 8):\n")
print(gene_tests_gse38129)

write.csv(gene_tests_gse38129,
          "GSE38129_gene_level_expression_tests_BH_adjusted.csv",
          row.names = FALSE)

# --- Individual gene paired Wilcoxon tests (Table 9) ------------------------
paired_gene_results <- lapply(genes_of_interest, function(gene) {
  normal_values <- gse38129_gene_long[
    gse38129_gene_long$tissue_type == "Normal", c("patient_id", gene)]
  tumor_values  <- gse38129_gene_long[
    gse38129_gene_long$tissue_type == "Tumor",  c("patient_id", gene)]
  colnames(normal_values)[2] <- "Normal"
  colnames(tumor_values)[2]  <- "Tumor"
  paired_values <- merge(normal_values, tumor_values,
                          by = "patient_id", all = FALSE)
  paired_values$Tumor_minus_Normal <- paired_values$Tumor - paired_values$Normal
  test_result <- wilcox.test(
    paired_values$Tumor, paired_values$Normal,
    paired = TRUE, exact = TRUE
  )
  data.frame(
    Gene = gene, N_pairs = nrow(paired_values),
    Normal_mean = mean(paired_values$Normal),
    Tumor_mean  = mean(paired_values$Tumor),
    Mean_Tumor_minus_Normal  = mean(paired_values$Tumor_minus_Normal),
    Median_Tumor_minus_Normal = median(paired_values$Tumor_minus_Normal),
    Paired_Wilcoxon_p_value  = test_result$p.value,
    Pairs_concordant = sum(paired_values$Tumor_minus_Normal > 0),
    Direction_in_tumor = ifelse(
      median(paired_values$Tumor_minus_Normal) > 0, "Higher", "Lower"),
    stringsAsFactors = FALSE
  )
})
paired_gene_results_gse38129 <- do.call(rbind, paired_gene_results)
paired_gene_results_gse38129$BH_FDR <- p.adjust(
  paired_gene_results_gse38129$Paired_Wilcoxon_p_value, method = "BH"
)
cat("\nGSE38129 Paired Gene Tests (Table 9):\n")
# Expected paired Wilcoxon p-values (exact = TRUE, matches manuscript Methods):
# SIM2   p=3.05e-5, BH-FDR=3.05e-5 (6/30 pairs Tumor>Normal — down-regulated)
# RFC4   p=5.59e-9, BH-FDR=1.40e-8 (29/30 pairs concordant)
# COL1A1 p=9.31e-9, BH-FDR=1.55e-8 (29/30 pairs concordant)
# MMP1   p=1.86e-9, BH-FDR=9.31e-9 (30/30 pairs concordant)
# CST1   p=4.66e-8, BH-FDR=5.82e-8 (27/30 pairs concordant)
print(paired_gene_results_gse38129)

write.csv(paired_gene_results_gse38129,
          "GSE38129_paired_gene_tests_BH.csv", row.names = FALSE)

# --- Paired composite score analysis -----------------------------------------
pair_audit_gse38129 <- gse38129_merged[
  , c("patient_id", "tissue_group", "coef_fixed_score")]

paired_wide_gse38129 <- reshape(
  pair_audit_gse38129,
  idvar = "patient_id", timevar = "tissue_group", direction = "wide"
)
paired_wide_gse38129$Tumor_minus_Normal <-
  paired_wide_gse38129$coef_fixed_score.Tumor -
  paired_wide_gse38129$coef_fixed_score.Normal

cat("\nPaired score summary (GSE38129):\n")
print(summary(paired_wide_gse38129$Tumor_minus_Normal))
cat("Pairs with higher tumour score:",
    sum(paired_wide_gse38129$Tumor_minus_Normal > 0), "/ 30\n")

# Two-sided exact paired Wilcoxon (primary reported result; matches
# manuscript Methods, which states an exact test was used)
wilcox_score_gse38129 <- wilcox.test(
  paired_wide_gse38129$coef_fixed_score.Tumor,
  paired_wide_gse38129$coef_fixed_score.Normal,
  paired = TRUE, alternative = "two.sided", exact = TRUE
)
cat("Paired Wilcoxon p =",
    format(wilcox_score_gse38129$p.value, scientific = TRUE), "\n")
# Expected: p = 1.86e-9
# Note: with exact = FALSE (normal approximation with continuity
# correction), this test instead returns p = 1.83e-6 — a materially
# different value. The exact test is reported in the manuscript
# because Methods Section 2.5 specifies an exact test throughout.

# Mean difference with 95% CI
mean_diff  <- mean(paired_wide_gse38129$Tumor_minus_Normal)
sd_diff    <- sd(paired_wide_gse38129$Tumor_minus_Normal)
n_pairs    <- nrow(paired_wide_gse38129)
se_diff    <- sd_diff / sqrt(n_pairs)
t_crit     <- qt(0.975, df = n_pairs - 1)
ci_lower   <- mean_diff - t_crit * se_diff
ci_upper   <- mean_diff + t_crit * se_diff
cat("Mean paired difference:", round(mean_diff, 2),
    "(95% CI", round(ci_lower, 2), "-", round(ci_upper, 2), ")\n")
# Expected: mean=5.25, 95% CI 4.45-6.05

write.csv(gse38129_merged,
          "GSE38129_coefficient_fixed_transfer_data.csv", row.names = FALSE)

# --- Paired score plot (Fig 5) -----------------------------------------------
paired_plot_data <- paired_wide_gse38129[
  order(paired_wide_gse38129$patient_id), ]

png("Fig5_GSE38129_paired_scores.png",
    width = 2100, height = 1500, res = 300)
plot(
  x = c(1, 2),
  y = range(paired_plot_data$coef_fixed_score.Normal,
             paired_plot_data$coef_fixed_score.Tumor),
  type = "n", xaxt = "n", xlab = "",
  ylab = "Coefficient-fixed, cohort-standardized score",
  main = "GSE38129: Matched Tumour\u2013Normal Score Differences"
)
axis(1, at = c(1, 2), labels = c("Adjacent normal", "ESCC tumour"))
abline(h = 0, lty = 2, col = "grey75")
for (i in seq_len(nrow(paired_plot_data))) {
  segments(x0 = 1, y0 = paired_plot_data$coef_fixed_score.Normal[i],
           x1 = 2, y1 = paired_plot_data$coef_fixed_score.Tumor[i],
           col = "#9E9E9E", lwd = 1)
}
points(rep(1, nrow(paired_plot_data)),
       paired_plot_data$coef_fixed_score.Normal,
       pch = 16, col = "#377EB8", cex = 1)
points(rep(2, nrow(paired_plot_data)),
       paired_plot_data$coef_fixed_score.Tumor,
       pch = 16, col = "#E41A1C", cex = 1)
legend("topleft",
       legend = c("Adjacent normal", "ESCC tumour"),
       pch = c(16, 16), col = c("#377EB8", "#E41A1C"),
       bty = "n", cex = 0.9)
usr <- par("usr")
# Positioned at x=1.15 / 0.30-0.38 (rather than 1.05 / 0.08-0.16) so the
# annotation text clears the legend box above it instead of overlapping
text(x = 1.15, y = usr[4] - 0.30 * diff(usr[3:4]),
     labels = "30/30 tumour scores higher", adj = 0, cex = 0.9)
text(x = 1.15, y = usr[4] - 0.38 * diff(usr[3:4]),
     labels = paste0("Paired Wilcoxon p = ",
                      format(wilcox_score_gse38129$p.value,
                             digits = 3, scientific = TRUE)),
     adj = 0, cex = 0.9)
dev.off()
# Label now pulled directly from wilcox_score_gse38129 (computed above
# with exact = TRUE), so the figure can never drift out of sync with
# the printed/reported p-value again.

# =============================================================================
# 7. TCGA EXPLORATORY DIAGNOSTIC ANALYSIS
# =============================================================================
# Note: Only 3 solid-tissue normal samples available in TCGA ESCC
# This analysis is exploratory only

# Load diagnostic subset (loaded from your TCGA data processing)
# merged_data_diag contains both -01 (tumour) and -11 (normal) samples

# Within diagnostic cohort z-standardization
for (gene in genes_of_interest) {
  merged_data_diag[[paste0("z_", gene)]] <-
    as.numeric(scale(as.numeric(merged_data_diag[[gene]])))
}

merged_data_diag$coef_fixed_score <-
  (published_coefficients["SIM2"]   * merged_data_diag$z_SIM2)   +
  (published_coefficients["RFC4"]   * merged_data_diag$z_RFC4)   +
  (published_coefficients["COL1A1"] * merged_data_diag$z_COL1A1) +
  (published_coefficients["MMP1"]   * merged_data_diag$z_MMP1)   +
  (published_coefficients["CST1"]   * merged_data_diag$z_CST1)

# --- Safety check: confirm tissue_type was built from strict -01/-11
# barcode matching only (Tumor=95, Normal=3), not a looser filter that
# could let a non-normal sample (e.g. a -06 metastatic barcode) through
# as "Normal". This exact mistake produced an incorrect Fig 2 during
# manuscript QC, where a 4th sample inflated the normal group to n=4
# and visibly distorted the ROC curve's step pattern. Fix before
# proceeding if this does not print Tumor=95, Normal=3.
cat("\nTCGA diagnostic subset tissue_type check:\n")
print(table(merged_data_diag$tissue_type))
stopifnot(
  "merged_data_diag must have exactly 95 Tumor samples" =
    sum(merged_data_diag$tissue_type == "Tumor") == 95,
  "merged_data_diag must have exactly 3 Normal samples (not 4) - check for a stray non -01/-11 barcode such as -06 metastatic" =
    sum(merged_data_diag$tissue_type == "Normal") == 3
)

# --- Table 5: Mean expression by tissue type ---------------------------------
table5_means <- lapply(genes_of_interest, function(gene) {
  normal_mean <- mean(
    as.numeric(merged_data_diag[[gene]][
      merged_data_diag$tissue_type == "Normal"]), na.rm = TRUE)
  tumor_mean <- mean(
    as.numeric(merged_data_diag[[gene]][
      merged_data_diag$tissue_type == "Tumor"]), na.rm = TRUE)
  data.frame(Gene = gene,
             Normal_mean = round(normal_mean, 3),
             Tumor_mean  = round(tumor_mean, 3),
             stringsAsFactors = FALSE)
})
table5_means <- do.call(rbind, table5_means)
cat("\n=== Table 5: Mean expression by tissue type (TCGA) ===\n")
print(table5_means)
# Expected:
# SIM2:   Normal=2.037, Tumor=1.360
# RFC4:   Normal=1.940, Tumor=4.092
# COL1A1: Normal=5.377, Tumor=9.194
# MMP1:   Normal=0.763, Tumor=5.441
# CST1:   Normal=0.041, Tumor=4.871

write.csv(table5_means, "Table5_TCGA_mean_expression.csv", row.names = FALSE)

# --- Table 6: Welch t-tests TCGA tumour vs normal (exploratory) -------------
# Note: n=3 normal samples; results are descriptive only
table6_welch <- lapply(genes_of_interest, function(gene) {
  normal_vals <- as.numeric(merged_data_diag[[gene]][
    merged_data_diag$tissue_type == "Normal"])
  tumor_vals  <- as.numeric(merged_data_diag[[gene]][
    merged_data_diag$tissue_type == "Tumor"])
  test <- t.test(tumor_vals, normal_vals, var.equal = FALSE)
  data.frame(
    Gene      = gene,
    p_value   = test$p.value,
    Direction = ifelse(mean(tumor_vals) > mean(normal_vals),
                       "Higher in tumour", "Lower in tumour"),
    Significant = ifelse(test$p.value < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
})
table6_welch <- do.call(rbind, table6_welch)
cat("\n=== Table 6: Welch t-tests TCGA tumour vs normal (exploratory) ===\n")
print(table6_welch)
# Expected p-values:
# SIM2=0.5664 (NS), RFC4=0.0030, COL1A1=0.0988 (NS),
# MMP1=0.0047, CST1<2.2e-16

write.csv(table6_welch,
          "Table6_TCGA_Welch_tests_exploratory.csv", row.names = FALSE)

# --- Paired tumour-normal check (3 matched pairs) ----------------------------
# Verify 3 participant-matched pairs show higher tumour score
tcga_normal_samples <- merged_data_diag$sample[
  merged_data_diag$tissue_type == "Normal"]
tcga_tumor_matched <- sub("-11$", "-01", tcga_normal_samples)

cat("\nTCGA matched pairs check:\n")
cat("Normal samples:", tcga_normal_samples, "\n")
cat("Matched tumour:", tcga_tumor_matched, "\n")

pairs_in_data <- tcga_tumor_matched[
  tcga_tumor_matched %in% merged_data_diag$sample]
cat("Matched pairs found:", length(pairs_in_data), "/ 3\n")

# --- ROC analysis (Fig 2, exploratory) --------------------------------------
roc_obj <- roc(
  response  = merged_data_diag$tissue_type,
  predictor = merged_data_diag$coef_fixed_score,
  levels    = c("Normal", "Tumor"),
  ci = TRUE, ci.method = "delong"
)

set.seed(42)
ci_boot_tcga <- ci.auc(roc_obj, method = "bootstrap",
                        boot.n = 2000, stratified = TRUE)

cat("\nTCGA Exploratory AUC:", round(as.numeric(auc(roc_obj)), 4), "\n")
cat("DeLong 95% CI:", round(ci(roc_obj)[1], 4), "-",
    round(ci(roc_obj)[3], 4), "\n")
cat("Bootstrap 95% CI:", round(ci_boot_tcga[1], 4), "-",
    round(ci_boot_tcga[3], 4), "\n")
# Expected: AUC=0.9895, DeLong CI 0.9688-1.0000, Bootstrap CI 0.9684-1.0000

png("Fig2_TCGA_Exploratory_ROC.png", width = 1800, height = 1500, res = 300)
plot(roc_obj, col = "#1B7837", lwd = 3,
     legacy.axes = TRUE, xlim = c(1, 0), ylim = c(0, 1),
     main = "TCGA ESCC (Exploratory): Coefficient-Fixed Score")
abline(a = 0, b = 1, lty = 2, col = "grey70")
text(x = 0.55, y = 0.15, adj = 0, cex = 0.9,
     labels = paste0(
       "AUC = 0.9895\n95 tumour / 3 solid-tissue normal\n",
       "DeLong 95% CI: 0.9688-1.0000\n",
       "(Exploratory: n=3 normal)"
     ))
dev.off()

write.csv(merged_data_diag,
          "ESCC_verified_diagnostic_dataset_95_tumor_3_normal.csv",
          row.names = FALSE)

# =============================================================================
# 8. SESSION INFO
# =============================================================================

cat("\n=== Session Information ===\n")
print(sessionInfo())

sink("session_info.txt")
sessionInfo()
sink()

cat("\n=== Analysis Complete ===\n")
cat("All figures saved as PNG files.\n")
cat("All result tables saved as CSV files.\n")
