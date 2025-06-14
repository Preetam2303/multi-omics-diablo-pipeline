## Multi-Omics DIABLO Pipeline for Immune Evasion Biomarker Discovery
This repository contains a complete multi-omics integration pipeline for identifying predictive biomarkers of immune evasion using DIABLO from the mixOmics package. It supports transcriptomics, DNA methylation, and proteomics data, and incorporates machine learning and functional enrichment analysis (GSEA) to validate the findings.

Key Highlight: Combines DIABLO-based multi-omics integration with LASSO-selected features, trains interpretable ML models (XGBoost, RF, SVM), and validates them via GO/KEGG enrichment and SHAP-based interpretation.

Applications
Multi-omics integration for systems biology

Immune evasion biomarker discovery

Predictive model building using LASSO-selected omics features

Functional pathway validation (GO / KEGG / GSEA)

SHAP-based interpretation of model decisions

Repository Structure
File	Description
1_data_preprocessing.R	Transcriptomics QC, normalization (Seurat), and immune evasion scoring
2_diablo_model.R	DIABLO model setup, training, and visualization
3_diablo_tuning.R	Cross-validation and tuning of DIABLO model parameters
4_lasso_feature_selection.R	LASSO regression for feature selection from each omic
5_ml_modeling.R	Builds ML models (SVM, RF, XGBoost) using combined omics features
6_shap_analysis.R	SHAP-based model interpretation and top gene identification
7_pathway_enrichment.R	GO and KEGG pathway analysis for top SHAP genes

# Requirements
R ≥ 4.0

Key R packages:

Seurat

mixOmics

glmnet, caret, xgboost, randomForest

clusterProfiler, enrichplot, SHAPforxgboost

ggplot2, dplyr, data.table, etc.

## Install all dependencies:

R
Copy
Edit
install.packages(c("Seurat", "mixOmics", "glmnet", "caret", "xgboost", "randomForest", "clusterProfiler", "SHAPforxgboost"))
How to Run
Clone this repo:

bash
Copy
Edit
git clone https://github.com/Preetam2303/multi-omics-diablo-pipeline.git
cd multi-omics-diablo-pipeline
Launch R and run scripts in order:

R
Copy
Edit
source("1_data_preprocessing.R")
source("2_diablo_model.R")
source("3_diablo_tuning.R")
source("4_lasso_feature_selection.R")
source("5_ml_modeling.R")
source("6_shap_analysis.R")
source("7_pathway_enrichment.R")
## Sample Output Highlights
DIABLO component plots and loading scores

LASSO-selected gene lists from each block

ROC curves and AUC values from ML models

SHAP summary and dependence plots

KEGG/GO pathway barplots for top predictive genes

## Citation and Background
This project is inspired by integrative biomarker discovery approaches using DIABLO (Rohart et al., 2017) and interpretable ML models for omics-based prediction.

If you use this pipeline, please cite the following:

- Banerjee, P. (2025). *Multi-Omics Integration and SHAP-Based Biomarker Discovery in Immune Evasion*. [Preprint in preparation].
- Rohart, F., Gautier, B., Singh, A., & Lê Cao, K.-A. (2017). mixOmics: An R package for 'omics feature selection and multiple data integration. *PLoS Computational Biology*, 13(11), e1005752.
- Lundberg, S. M., & Lee, S.-I. (2017). A unified approach to interpreting model predictions. *Advances in Neural Information Processing Systems*, 30.
- Yu, G., Wang, L.-G., Han, Y., & He, Q.-Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS: A Journal of Integrative Biology*, 16(5), 284–287.

