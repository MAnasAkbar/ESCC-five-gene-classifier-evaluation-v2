# ESCC Five-Gene Classifier: Independent Evaluation

## Study

**Prognostic Evaluation and Multi-Cohort Cross-Platform Diagnostic Replication of a Five-Gene Transcriptomic Classifier in Esophageal Squamous Cell Carcinoma**

Muhammad Anas Akbar, Nawera Khan
Department of Biotechnology, Government College University, Lahore, Pakistan
Corresponding author: rnanasje@gmail.com

Prepared for submission to PLOS ONE, 2026.

---

## Version History

**v1.2 (current)** — Two further fixes found during figure QC against the
manuscript: (1) Fig 5's p-value annotation overlapped the plot legend;
repositioned so both are legible. (2) Added a `stopifnot()` safety check
before Table 5/Table 6 confirming the TCGA diagnostic subset has exactly
95 Tumor / 3 Normal samples. This guards against a real failure mode
found during QC, where a separate (non-script) interactive object had
let a `-06` (metastatic) barcode slip through as "Normal," inflating the
group to n=4 and visibly distorting that cohort's ROC curve. The
script's own filtering logic was already correct and unaffected; this
check only makes that guarantee explicit and loud if it's ever violated.

**v1.1** — Corrected the paired Wilcoxon signed-rank tests in the
GSE38129 analysis (individual genes, Table 9, and the composite score) to
use `exact = TRUE`, matching the exact two-sided test specified in the
manuscript's Methods (Section 2.5). The previous release used
`exact = FALSE` (a normal approximation with continuity correction),
which gave numerically different p-values than the stated method. See
`results/GSE38129_paired_gene_tests_BH.csv` and
`results/GSE38129_composite_score_paired_Wilcoxon.csv` for both the
corrected values and, for transparency, the original approximate values.
No other results, cohorts, or conclusions in the manuscript were
affected by this correction — GSE38129 remained significant for all
five genes and the composite score before and after the fix; only the
reported p-value magnitudes and the exact composite-score value changed
(from the mislabelled 9.31 × 10⁻¹⁰ to the correct exact two-sided value,
1.86 × 10⁻⁹).

**v1.0** — Initial deposit accompanying manuscript submission.

---

## Overview

This repository contains the R analysis code for an independent
evaluation of the five-gene ESCC transcriptomic classifier (SIM2, RFC4,
COL1A1, MMP1, CST1) originally published by Khalil et al. (Bioscience
Insights, 2026;3(2):478-495).

The published coefficients were applied unchanged (coefficient-fixed)
with within-cohort z-standardization across three publicly available
evaluation cohorts: TCGA-ESCA, GSE20347, and GSE38129.

---

## Data Sources

All data are publicly available. No new human subjects data were
generated.

| Cohort | Source | Accession | Access |
|--------|--------|-----------|--------|
| TCGA-ESCA | UCSC Xena | TCGA-ESCA | https://xenabrowser.net |
| GSE20347 | GEO | GSE20347 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20347 |
| GSE38129 | GEO | GSE38129 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38129 |

---

## Repository Structure

```
├── ESCC_analysis_complete.R    # Complete, single-file analysis script
│                                 (data acquisition through figures)
├── MANIFEST.txt                # GDC file manifest for TCGA-ESCA raw data
├── README.md
└── results/
    ├── ESCC_verified_diagnostic_dataset_95_tumor_3_normal.csv
    ├── GSE20347_verified_pair_mapping.csv
    ├── GSE20347_coefficient_fixed_transfer_data.csv
    ├── GSE20347_five_gene_paired_Wilcoxon_results.csv
    ├── GSE20347_signature_score_paired_differences.csv
    ├── GSE38129_coefficient_fixed_transfer_data.csv
    ├── GSE38129_gene_level_expression_tests.csv
    ├── GSE38129_gene_level_expression_tests_BH_adjusted.csv
    ├── GSE38129_paired_gene_tests_BH.csv
    ├── GSE38129_composite_score_paired_Wilcoxon.csv
    └── GSE121931_sample_ID_mapping.csv
```

`ESCC_analysis_complete.R` runs end to end in eight numbered sections
(setup, TCGA data and survival cohort, survival analysis, GEO data
acquisition, score calculation, GSE20347 analysis, GSE38129 analysis,
TCGA exploratory diagnostic analysis, session info). Each section
prints its key results to the console with the expected value given
as a comment, so output can be checked against the manuscript as the
script runs.

---

## How to Run

```r
source("ESCC_analysis_complete.R")
```

**Note:** TCGA-ESCA data must be downloaded manually from UCSC Xena
before the TCGA sections will run (the GDC manifest for the raw
RNA-seq files is provided in `MANIFEST.txt`). The GEO sections
(GSE20347, GSE38129) download automatically via `GEOquery` and need no
manual steps.

---

## R Environment

```
R version: 4.6.1
Packages: survival, survminer, pROC, GEOquery
```

Key package versions:
- survival 3.8
- survminer (2024)
- pROC 1.18
- GEOquery 2.70

---

## Published Coefficients Used

From Khalil et al. (2026), Table 3:

| Gene   | Coefficient |
|--------|-------------|
| SIM2   | -0.824      |
| RFC4   | +1.079      |
| COL1A1 | +0.639      |
| MMP1   | +0.712      |
| CST1   | +0.521      |

---

## Key Results (current, v1.1)

| Cohort | n (T/N) | AUC | Paired Wilcoxon p (composite score, exact two-sided) |
|--------|---------|-----|-------------------------------------------------------|
| TCGA (exploratory) | 95/3 | 0.9895 | — |
| GSE20347 | 17/17 | 1.000 | 1.53 × 10⁻⁵ |
| GSE38129 | 30/30 | 0.9611 | 1.86 × 10⁻⁹ |

TCGA prognostic: Cox HR = 0.962 (95% CI 0.817-1.132), p = 0.641
(n = 95, 32 events)

GSE38129 individual-gene paired Wilcoxon (exact two-sided; Table 9):

| Gene | Pairs concordant | p (exact) | BH-FDR |
|------|-------------------|-----------|--------|
| SIM2 | 6/30 | 3.05 × 10⁻⁵ | 3.05 × 10⁻⁵ |
| RFC4 | 29/30 | 5.59 × 10⁻⁹ | 1.40 × 10⁻⁸ |
| COL1A1 | 29/30 | 9.31 × 10⁻⁹ | 1.55 × 10⁻⁸ |
| MMP1 | 30/30 | 1.86 × 10⁻⁹ | 9.31 × 10⁻⁹ |
| CST1 | 27/30 | 4.66 × 10⁻⁸ | 5.82 × 10⁻⁸ |

---

## License

Code: MIT License
Data: Subject to original data source terms (TCGA, GEO - publicly
available)

---

## Citation

If you use this code, please cite:

> Akbar MA, Khan N. Prognostic Evaluation and Multi-Cohort
> Cross-Platform Diagnostic Replication of a Five-Gene Transcriptomic
> Classifier in Esophageal Squamous Cell Carcinoma. PLOS ONE. 2026.
> [DOI pending]

---

## Results Files Description

| File | Description |
|------|-------------|
| ESCC_verified_diagnostic_dataset_95_tumor_3_normal.csv | TCGA ESCC diagnostic subset (95 tumour, 3 normal) |
| GSE20347_verified_pair_mapping.csv | Confirmed tumour-normal pair IDs for GSE20347 |
| GSE20347_coefficient_fixed_transfer_data.csv | Per-sample scores and z-scores for GSE20347 |
| GSE20347_five_gene_paired_Wilcoxon_results.csv | Individual gene paired Wilcoxon results (GSE20347) |
| GSE20347_signature_score_paired_differences.csv | Paired score differences (Normal vs Tumour), GSE20347 |
| GSE38129_coefficient_fixed_transfer_data.csv | Per-sample scores and z-scores for GSE38129 |
| GSE38129_gene_level_expression_tests.csv | Initial Welch t-test results, GSE38129 |
| GSE38129_gene_level_expression_tests_BH_adjusted.csv | Welch tests with BH-FDR correction (Table 8) |
| GSE38129_paired_gene_tests_BH.csv | Exact paired Wilcoxon results with BH-FDR (Table 9; corrected in v1.1) |
| GSE38129_composite_score_paired_Wilcoxon.csv | Composite-score paired Wilcoxon test, exact vs. approximate (added in v1.1) |
| GSE121931_sample_ID_mapping.csv | Sample ID mapping file |

---

## Notes on MANIFEST.txt

The MANIFEST.txt file contains the GDC (Genomic Data Commons) file
manifest used to download TCGA-ESCA RNA-seq files. These files are not
included in the repository due to size but can be re-downloaded from
GDC using this manifest at: https://portal.gdc.cancer.gov/
