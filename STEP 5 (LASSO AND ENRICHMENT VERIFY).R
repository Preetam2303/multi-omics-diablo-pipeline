load("diablo_pipeline_results1.RData")
load("diablo_pipeline_results1_with_feats.RData")
load("selected_features.RData")
####################
print(prot_feats_loadings)
ls()
#######################
library(glmnet)
library(caret)

# Assume these exist from your DIABLO step:
# prot_feats_loadings, meth_feats_loadings, rna_feats_loadings
# X_train, Y_train, X_test, Y_test

# Prepare a named list of features from loadings
selected_features <- list(
  proteomics = prot_feats_loadings,
  methylation = meth_feats_loadings,
  transcriptomics = rna_feats_loadings
)

# Initialize lists to store results
lasso_models <- list()
test_preds <- list()
test_probs <- list()
test_cm <- list()

set.seed(42)

for (block in names(selected_features)) {
  
  cat(sprintf("Running LASSO for block: %s\n", block))
  
  # Subset training data to selected features
  feats <- selected_features[[block]]
  
  # Make sure features exist in the data (intersect)
  feats <- intersect(feats, colnames(X_train[[block]]))
  
  X_tr <- X_train[[block]][, feats, drop = FALSE]
  X_te <- X_test[[block]][, feats, drop = FALSE]
  
  # glmnet expects matrix input
  X_tr_mat <- as.matrix(X_tr)
  X_te_mat <- as.matrix(X_te)
  
  # Ensure Y_train is a factor
  y_tr <- factor(Y_train)
  
  # Convert factor to binary numeric response for glmnet (E.g., Evasive=1, Non-Evasive=0)
  y_tr_bin <- ifelse(y_tr == "Evasive", 1, 0)
  
  # Fit LASSO logistic regression with 10-fold CV to select lambda
  cvfit <- cv.glmnet(
    x = X_tr_mat,
    y = y_tr_bin,
    family = "binomial",
    alpha = 1,
    nfolds = 10,
    type.measure = "class"
  )
  
  best_lambda <- cvfit$lambda.min
  cat(sprintf("Best lambda for %s: %f\n", block, best_lambda))
  
  # Fit final model at best lambda
  lasso_model <- glmnet(
    x = X_tr_mat,
    y = y_tr_bin,
    family = "binomial",
    alpha = 1,
    lambda = best_lambda
  )
  
  lasso_models[[block]] <- lasso_model
  
  # Predict on test set (probabilities)
  probs <- predict(lasso_model, newx = X_te_mat, type = "response")
  test_probs[[block]] <- probs
  
  # Predict class labels (threshold 0.5)
  preds <- ifelse(probs >= 0.5, "Evasive", "Non-Evasive")
  preds <- factor(preds, levels = levels(y_tr))
  
  test_preds[[block]] <- preds
  
  # Confusion matrix
  cm <- confusionMatrix(preds, factor(Y_test, levels = levels(y_tr)))
  test_cm[[block]] <- cm
  
  cat(sprintf("\nConfusion Matrix for LASSO on %s:\n", block))
  print(cm)
  cat("\n----------------------------------------\n")
}

# Now lasso_models, test_preds, test_probs, and test_cm hold results per block
save(lasso_models, test_preds, test_probs, test_cm, file = "lasso_models_results.RData")

selected_lasso_features <- list()

for (block in names(lasso_models)) {
  coef_mat <- coef(lasso_models[[block]])
  # Remove intercept and get features with non-zero coefficients
  selected_features <- rownames(coef_mat)[which(coef_mat[,1] != 0)]
  selected_features <- setdiff(selected_features, "(Intercept)")
  
  selected_lasso_features[[block]] <- selected_features
  cat(sprintf("LASSO selected features for %s (%d):\n", block, length(selected_features)))
  print(selected_features)
}
save(selected_lasso_features, file = "lasso_selected_features.RData")
##################
library(clusterProfiler)
library(org.Hs.eg.db)

# Example: Convert gene symbols to Entrez IDs
convert_to_entrez <- function(gene_symbols) {
  bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}
# For proteomics block
prot_genes <- selected_lasso_features$proteomics
prot_entrez <- convert_to_entrez(prot_genes)$ENTREZID

# Perform GO Biological Process enrichment
prot_go <- enrichGO(gene         = prot_entrez,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    pAdjustMethod= "BH",
                    qvalueCutoff = 0.05,
                    readable     = TRUE)

# View top pathways
head(prot_go)
#
library(enrichplot)

dotplot(prot_go, showCategory = 10) + ggtitle("Proteomics GO BP Enrichment")
######################################
# Original proteomics genes
prot_genes <- selected_lasso_features$proteomics

# Convert to Entrez IDs
prot_entrez_df <- convert_to_entrez(prot_genes)

# Check mapping success
cat("Number of original proteomics genes:", length(prot_genes), "\n")
cat("Number of genes successfully mapped to Entrez IDs:", nrow(prot_entrez_df), "\n")

# If you want to see which genes failed mapping:
unmapped_prot_genes <- setdiff(prot_genes, prot_entrez_df$SYMBOL)
if(length(unmapped_prot_genes) > 0) {
  cat("Unmapped proteomics genes:\n")
  print(unmapped_prot_genes)
} else {
  cat("All proteomics genes were mapped successfully.\n")
}

######################################################
# Extract transcriptomics genes
transcript_genes <- selected_lasso_features$transcriptomics

# Convert transcriptomics genes to Entrez IDs
transcript_entrez_df <- convert_to_entrez(transcript_genes)

# Check mapping success
cat("Number of original transcriptomics genes:", length(transcript_genes), "\n")
cat("Number of genes successfully mapped to Entrez IDs:", nrow(transcript_entrez_df), "\n")

# List unmapped genes if any
unmapped_transcript_genes <- setdiff(transcript_genes, transcript_entrez_df$SYMBOL)
if(length(unmapped_transcript_genes) > 0) {
  cat("Unmapped transcriptomics genes:\n")
  print(unmapped_transcript_genes)
} else {
  cat("All transcriptomics genes were mapped successfully.\n")
}
###################
library(org.Hs.eg.db)
ensembl_id <- "ENSG00000006282"

symbol_df <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = ensembl_id,
                                   columns = c("SYMBOL", "ENTREZID"),
                                   keytype = "ENSEMBL")

print(symbol_df)
#############
# Replace Ensembl ID with symbol
transcript_genes <- gsub("ENSG00000006282", "SPATA20", transcript_genes)

# Convert gene symbols to Entrez IDs again
transcript_entrez_df <- bitr(transcript_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform GO enrichment as before
if (!is.null(transcript_entrez_df) && nrow(transcript_entrez_df) > 0) {
  
  enrich_res <- enrichGO(
    gene          = transcript_entrez_df$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  print(head(enrich_res))
  
  dotplot(enrich_res, showCategory = 10) + ggtitle("Transcriptomics GO BP Enrichment")
  
} else {
  cat("No valid gene IDs found for enrichment in transcriptomics block.\n")
}
print(enrich_res)
#run this 
library(clusterProfiler)
kegg_res <- enrichKEGG(
  gene         = transcript_entrez_df$ENTREZID,
  organism     = 'hsa',
  pvalueCutoff = 0.2
)
print(kegg_res)
dotplot(kegg_res, showCategory = 10) + ggtitle("Transcriptomics KEGG Pathway Enrichment")
#######################
####################################
# Methylation block KEGG enrichment
####################################

# Make sure the following libraries are loaded
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Extract methylation block genes
methy_genes <- selected_lasso_features[["methylation"]]

# Convert methylation gene symbols to Entrez IDs
methy_entrez_df <- tryCatch({
  bitr(methy_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, error = function(e) {
  message("Conversion failed for methylation block.")
  NULL
})

# Check mapping success
if (!is.null(methy_entrez_df) && nrow(methy_entrez_df) > 0) {
  cat("Number of original methylation genes:", length(methy_genes), "\n")
  cat("Number of genes successfully mapped to Entrez IDs:", nrow(methy_entrez_df), "\n")
  
  unmapped_methy_genes <- setdiff(methy_genes, methy_entrez_df$SYMBOL)
  if (length(unmapped_methy_genes) > 0) {
    cat("Unmapped methylation genes:\n")
    print(unmapped_methy_genes)
  } else {
    cat("All methylation genes were mapped successfully.\n")
  }
  
  # Perform KEGG enrichment
  kegg_methy_res <- enrichKEGG(
    gene         = methy_entrez_df$ENTREZID,
    organism     = 'hsa',
    pvalueCutoff = 0.2
  )
  
  # Print and plot results
  print(kegg_methy_res)
  
  if (nrow(kegg_methy_res) > 0) {
    dotplot(kegg_methy_res, showCategory = 10) + ggtitle("Methylation Block KEGG Pathway Enrichment")
  } else {
    cat("No significant KEGG pathways found for methylation block.\n")
  }
  
} else {
  cat("No valid gene IDs found for KEGG enrichment in methylation block.\n")
}
# Save proteomics GO enrichment
save(prot_go, prot_entrez_df, file = "proteomics_go_enrichment.RData")

# Save transcriptomics GO and KEGG enrichment
save(enrich_res, kegg_res, transcript_entrez_df, file = "transcriptomics_enrichment.RData")

list.files(pattern = "enrichment.RData")

