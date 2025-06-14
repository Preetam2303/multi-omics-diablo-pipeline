# Load the required data objects
load("lasso_models_results.RData")
load("lasso_selected_features.RData")
load("diablo_pipeline_results1_with_feats.RData")
load("selected_features.RData")
# Combine features from all omics blocks
all_selected_features <- unique(unlist(selected_lasso_features))
length(all_selected_features)
print(all_selected_features)
###############
str(X_train)
#
# Extract each block
X_train_prot <- X_train$proteomics
X_train_meth <- X_train$methylation
X_train_rna  <- X_train$transcriptomics

X_test_prot <- X_test$proteomics
X_test_meth <- X_test$methylation
X_test_rna  <- X_test$transcriptomics

# Subset based on selected features from each block
X_train_sub <- cbind(
  X_train_prot[, intersect(colnames(X_train_prot), all_selected_features), drop=FALSE],
  X_train_meth[, intersect(colnames(X_train_meth), all_selected_features), drop=FALSE],
  X_train_rna[, intersect(colnames(X_train_rna), all_selected_features), drop=FALSE]
)

X_test_sub <- cbind(
  X_test_prot[, intersect(colnames(X_test_prot), all_selected_features), drop=FALSE],
  X_test_meth[, intersect(colnames(X_test_meth), all_selected_features), drop=FALSE],
  X_test_rna[, intersect(colnames(X_test_rna), all_selected_features), drop=FALSE]
)

# Confirm dimensions
dim(X_train_sub)
dim(X_test_sub)
#################
# Convert Y vectors to factor
Y_train_factor <- as.factor(Y_train)
Y_test_factor  <- as.factor(Y_test)
####################
library(randomForest)

# Train
set.seed(123)
rf_model <- randomForest(x = X_train_sub, y = Y_train_factor, ntree = 500, importance = TRUE)

# Predict
rf_preds <- predict(rf_model, X_test_sub)
rf_probs <- predict(rf_model, X_test_sub, type = "prob")

# Confusion Matrix
rf_cm <- table(Predicted = rf_preds, Actual = Y_test_factor)
print(rf_cm)
##############################
library(xgboost)

# Convert to matrix for xgboost
dtrain <- xgb.DMatrix(data = as.matrix(X_train_sub), label = as.numeric(Y_train_factor) - 1)
dtest  <- xgb.DMatrix(data = as.matrix(X_test_sub),  label = as.numeric(Y_test_factor) - 1)

# Train
set.seed(123)
xgb_model <- xgboost(data = dtrain, nrounds = 100, objective = "binary:logistic", verbose = 0)

# Predict
xgb_probs <- predict(xgb_model, dtest)
xgb_preds <- ifelse(xgb_probs > 0.5, 1, 0)

# Confusion Matrix
xgb_cm <- table(Predicted = xgb_preds, Actual = as.numeric(Y_test_factor) - 1)
print(xgb_cm)
###############################
# Accuracy
rf_acc <- mean(rf_preds == Y_test_factor)
xgb_acc <- mean(xgb_preds == (as.numeric(Y_test_factor) - 1))
cat("Random Forest Accuracy:", rf_acc, "\n")
cat("XGBoost Accuracy:", xgb_acc, "\n")

# AUC
library(pROC)
rf_auc <- roc(Y_test_factor, rf_probs[, 2])
xgb_auc <- roc(Y_test_factor, xgb_probs)
cat("RF AUC:", auc(rf_auc), "\n")
cat("XGB AUC:", auc(xgb_auc), "\n")
##########
library(e1071)
set.seed(123)
svm_model <- svm(x = X_train_sub, y = Y_train_factor, 
                 kernel = "radial", probability = TRUE, 
                 scale = TRUE)
svm_preds <- predict(svm_model, X_test_sub)
svm_probs <- attr(predict(svm_model, X_test_sub, probability = TRUE), "probabilities")[,2]
##
svm_cm <- table(Predicted = svm_preds, Actual = Y_test_factor)
print(svm_cm)

svm_acc <- mean(svm_preds == Y_test_factor)
cat("SVM Accuracy:", svm_acc, "\n")
#
library(pROC)
svm_auc <- roc(Y_test_factor, svm_probs)
cat("SVM AUC:", auc(svm_auc), "\n")
#################
install.packages(c("pROC", "ggplot2", "caret", "yardstick", "reshape2"))
library(pROC)
library(ggplot2)
library(caret)
library(yardstick)
library(reshape2)
##############
# Create ROC curves
rf_roc  <- roc(Y_test_factor, rf_probs[, 2])
xgb_roc <- roc(Y_test_factor, xgb_probs)
svm_roc <- roc(Y_test_factor, svm_probs)

# Plot
plot(rf_roc, col = "steelblue", lwd = 2, main = "ROC Curves", legacy.axes = TRUE)
plot(xgb_roc, col = "darkgreen", lwd = 2, add = TRUE)
plot(svm_roc, col = "darkred", lwd = 2, add = TRUE)
legend("bottomright", legend = c("Random Forest", "XGBoost", "SVM"),
       col = c("steelblue", "darkgreen", "darkred"), lwd = 2)
################
# Load necessary libraries
library(caret)
library(ggplot2)
library(reshape2)

# Function to compute evaluation metrics
compute_metrics <- function(preds, probs, actual) {
  # Ensure predictions have the same factor levels as actual
  preds <- factor(preds, levels = levels(actual))
  
  cm <- confusionMatrix(preds, actual)
  
  acc <- cm$overall["Accuracy"]
  f1 <- ifelse(is.na(cm$byClass["F1"]), 0, cm$byClass["F1"])
  prec <- ifelse(is.na(cm$byClass["Precision"]), 0, cm$byClass["Precision"])
  rec <- ifelse(is.na(cm$byClass["Recall"]), 0, cm$byClass["Recall"])
  rmsd <- sqrt(mean((as.numeric(actual) - as.numeric(preds))^2))
  
  return(c(Accuracy = acc, F1 = f1, Precision = prec, Recall = rec, RMSD = rmsd))
}

# Ensure predicted labels are factors with correct levels
xgb_preds_factor <- factor(xgb_preds, levels = c(0, 1), labels = levels(Y_test_factor))
rf_preds_factor  <- factor(rf_preds,  levels = levels(Y_test_factor))  # Optional, just to be safe

# Recalculate metrics
rf_metrics  <- compute_metrics(rf_preds_factor, rf_probs[, 2], Y_test_factor)
xgb_metrics <- compute_metrics(xgb_preds_factor, xgb_probs, Y_test_factor)
svm_metrics <- compute_metrics(svm_preds, svm_probs, Y_test_factor)

# Combine metrics into a clean dataframe
metrics_df <- data.frame(
  Metric = c("Accuracy", "F1 Score", "Precision", "Recall", "RMSD"),
  RandomForest = as.numeric(rf_metrics),
  XGBoost = as.numeric(xgb_metrics),
  SVM = as.numeric(svm_metrics)
)

# Melt for plotting
metrics_melt <- melt(metrics_df, id.vars = "Metric", variable.name = "Model", value.name = "Value")

# Plot: Publication-quality performance bar plot
ggplot(metrics_melt, aes(x = Metric, y = Value, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  scale_fill_manual(values = c("steelblue", "darkgreen", "darkred")) +
  theme_minimal(base_size = 14) +
  labs(title = "Model Performance Comparison",
       y = "Metric Value",
       x = "Metric",
       fill = "Model") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

# Print SVM metrics separately if needed
print(svm_metrics)
#
levels(Y_test_factor)
levels(svm_preds)
levels(rf_preds)
levels(xgb_preds)  # If it's numeric, then convert
###############
print(xgb_metrics)
print(rf_metrics)
#######################################SHAP START
# Install if not already installed
install.packages("iml")
install.packages("mlr")     # Required for creating Predictor object
install.packages("data.table")

# Load libraries
library(iml)
library(mlr)
library(data.table)
library(xgboost)
#
# Convert to matrix and assign column names BEFORE creating DMatrix
train_matrix <- as.matrix(X_train_sub)
colnames(train_matrix) <- colnames(X_train_sub)  # Set feature names here

# Now create DMatrix with column names already set
dtrain <- xgb.DMatrix(data = train_matrix, label = as.numeric(Y_train_factor) - 1)

# Train XGBoost model
xgb_model <- xgboost(
  data = dtrain,
  objective = "binary:logistic",
  nrounds = 100,
  verbose = 0
)

# Align test set to training feature names
X_test_df <- as.data.frame(X_test_sub)
X_test_df <- X_test_df[, colnames(X_train_sub), drop = FALSE]  # Ensure same columns and order

########
predict_function <- function(model, newdata) {
  test_matrix <- as.matrix(newdata)
  colnames(test_matrix) <- colnames(X_train_sub)  # Align feature names here too
  dnew <- xgb.DMatrix(data = test_matrix)
  predict(model, dnew)
}
########
ncol(X_test_df) == length(colnames(X_train_sub))  # Should return TRUE
predict_function <- function(model, newdata) {
  test_matrix <- as.matrix(newdata)
  
  # Check and fix only if dimensions match
  if (ncol(test_matrix) == length(colnames(X_train_sub))) {
    colnames(test_matrix) <- colnames(X_train_sub)
  }
  
  dnew <- xgb.DMatrix(data = test_matrix)
  predict(model, dnew)
}
###
predictor <- Predictor$new(
  model = xgb_model,
  data = X_test_df,
  y = Y_test_factor,
  predict.function = predict_function,
  type = "prob"
)

feature_imp <- FeatureImp$new(predictor, loss = "ce")
plot(feature_imp)
##

print(colnames(X_train_sub))
print(colnames(X_test_df))
################
# Ensure unique and matching column names in both training and test sets
colnames(X_train_sub) <- make.names(colnames(X_train_sub), unique = TRUE)
colnames(X_test_sub)  <- make.names(colnames(X_test_sub), unique = TRUE)

# Convert them to data frames if needed
X_test_df <- as.data.frame(X_test_sub)
X_test_df <- X_test_df[, colnames(X_train_sub), drop = FALSE]

#
dtrain <- xgb.DMatrix(data = as.matrix(X_train_sub), label = as.numeric(Y_train_factor) - 1)
dtest  <- xgb.DMatrix(data = as.matrix(X_test_df), label = as.numeric(Y_test_factor) - 1)

set.seed(123)
xgb_model <- xgboost(data = dtrain, nrounds = 100, objective = "binary:logistic", verbose = 0)


#
predictor <- Predictor$new(
  model = xgb_model,
  data = X_test_df,
  y = Y_test_factor,
  predict.function = predict_function,
  type = "prob"
)
feature_imp <- FeatureImp$new(predictor, loss = "ce")
plot(feature_imp)
###########
# Compute mean absolute SHAP values for each feature
shap_means <- colMeans(abs(shap_values))

# Create a data frame for plotting
shap_df <- data.frame(
  Feature = names(shap_means),
  MeanAbsSHAP = shap_means
)

# Load ggplot2 (if not already loaded)
library(ggplot2)

# Plot: a horizontal bar chart sorted by importance
ggplot(shap_df, aes(x = reorder(Feature, MeanAbsSHAP), y = MeanAbsSHAP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Feature", y = "Mean |SHAP| Value",
       title = "SHAP Value Summary")
#############################
# Compute mean absolute SHAP values for each feature
shap_means <- colMeans(abs(shap_values))

# Select the top 15 features based on the mean absolute SHAP values
top_features <- sort(shap_means, decreasing = TRUE)[1:20]

# Create a data frame for these top features
shap_df <- data.frame(
  Feature = names(top_features),
  MeanAbsSHAP = as.numeric(top_features)
)

# Load ggplot2 if not already loaded
library(ggplot2)

# Plot: a horizontal bar chart sorted by importance
ggplot(shap_df, aes(x = reorder(Feature, MeanAbsSHAP), y = MeanAbsSHAP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Feature", y = "Mean |SHAP| Value",
       title = "Top 15 SHAP Value Summary")

#############################################
library(ggplot2)

# Identify top features based on mean absolute SHAP values
shap_means <- colMeans(abs(shap_values))
top_features <- names(sort(shap_means, decreasing = TRUE))[5]  # top 3 features

# Loop over the top features and plot dependence plots
for (feat in top_features) {
  # Create a data frame with the feature values and its corresponding SHAP values
  df <- data.frame(
    FeatureValue = as.numeric(X_test_df[, feat]),
    SHAP = shap_values[, feat]
  )
  
  p <- ggplot(df, aes(x = FeatureValue, y = SHAP)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_smooth(method = "loess", se = FALSE, color = "darkred") +
    labs(
      title = paste("Dependence Plot for", feat),
      x = feat,
      y = "SHAP value"
    ) +
    theme_minimal(base_size = 14)
  
  print(p)
}

####################
# Install SHAPforxgboost if it's not already installed
if (!requireNamespace("SHAPforxgboost", quietly = TRUE)) {
  install.packages("SHAPforxgboost")
}
# Load the package
library(SHAPforxgboost)

# Compute SHAP values using your trained XGBoost model and training data.
# We use as.matrix() to convert your training data into matrix format.
shap_values_xgb <- shap.values(xgb_model = xgb_model, X_train = as.matrix(X_train_sub))

# Prepare the long-format data required for plotting.
shap_long <- shap.prep(shap_contrib = shap_values_xgb$shap_score, 
                       X_train = as.matrix(X_train_sub))

# Create the force plot for a selected observation (for example, row_index = 1)
# You can change 'row_index' to view other observations.
SHAPforxgboost:::shap.plot.force(data_long = shap_long, row_index = 1)
#############
if (!requireNamespace("shapviz", quietly = TRUE)) {
  install.packages("shapviz")
}
library(shapviz)
sv_obj <- shapviz(xgb_model,
                  X_train = as.matrix(X_train_sub),
                  X_pred  = as.matrix(X_test_df))

sv_force(sv_obj, row = 1)
######
# --------------------------------------------------------------------------------
# Step 1: Identify Top 15 Features by SHAP Importance
# --------------------------------------------------------------------------------
# Assume 'shap_values' is your matrix of SHAP values with gene names as column names.
# Calculate mean absolute SHAP values and select the top 15 features.
shap_means <- colMeans(abs(shap_values))
top15_genes <- names(sort(shap_means, decreasing = TRUE))[1:15]
cat("Top 15 genes:\n")
print(top15_genes)

# --------------------------------------------------------------------------------
# Step 2: Map Gene Symbols to Entrez IDs
# --------------------------------------------------------------------------------
# Install and load required packages if not already installed
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  install.packages("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  install.packages("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)

# Convert the top 15 gene symbols to Entrez IDs.
top15_mapping <- bitr(top15_genes,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
cat("Mapping of top 15 genes:\n")
print(top15_mapping)

# --------------------------------------------------------------------------------
# Step 3: Perform Enrichment Analysis (KEGG Over-Representation)
# --------------------------------------------------------------------------------
# Using the enrichKEGG function to look for KEGG pathways enriched in the Top 15 genes.
enrich_kegg <- enrichKEGG(gene         = top15_mapping$ENTREZID,
                          organism     = "hsa",       # 'hsa' for human
                          pAdjustMethod = "BH",
                          qvalueCutoff  = 0.05)

# View the enrichment results as a data frame
enrich_results <- as.data.frame(enrich_kegg)
cat("Enrichment results (top rows):\n")
head(enrich_results)

# --------------------------------------------------------------------------------
# Optional: Visualize the Enrichment Results
# --------------------------------------------------------------------------------
# A dot plot can nicely summarize the enriched pathways.
dotplot(enrich_kegg, showCategory = 15, title = "KEGG Pathway Enrichment of Top 15 Genes")
############################################
# --------------------------------------------------------------------------------
# Step 1: Create a Ranked Gene List from SHAP Values
# --------------------------------------------------------------------------------
shap_means <- colMeans(abs(shap_values))
# Rank all genes in decreasing order:
gene_rank <- sort(shap_means, decreasing = TRUE)
cat("Total genes in ranked list:", length(gene_rank), "\n")

# --------------------------------------------------------------------------------
# Step 2: Map Gene Symbols to Entrez IDs for the Full List  
# --------------------------------------------------------------------------------
all_mapping <- bitr(names(gene_rank),
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
# Keep only the genes that were successfully mapped
mapped_genes <- all_mapping$SYMBOL
gene_rank <- gene_rank[names(gene_rank) %in% mapped_genes]

# Replace gene symbols with their corresponding ENTREZ IDs
names(gene_rank) <- all_mapping$ENTREZID[match(names(gene_rank), all_mapping$SYMBOL)]
cat("Ranked gene list prepared for GO analysis. Length:", length(gene_rank), "\n")
# Now, gene_rank should contain the 27 mapped genes.

# --------------------------------------------------------------------------------
# Step 3: Run GO Enrichment Analysis (ORA) Using enrichGO()
# --------------------------------------------------------------------------------
# In this example, we'll analyze Biological Process (BP).
ego <- enrichGO(gene         = names(gene_rank),
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",   # our gene_rank names are ENTREZ IDs now
                ont          = "BP",        # Change to "MF" or "CC" for Molecular Function or Cellular Component, if desired.
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)        # converts Entrez IDs to Gene Symbols for easier interpretation

# Convert the enrichment results to a data frame and show the top results.
ego_results <- as.data.frame(ego)
cat("GO Enrichment results (top rows):\n")
head(ego_results)

# --------------------------------------------------------------------------------
# Optional: Visualize the GO Enrichment Results
# --------------------------------------------------------------------------------
# Generate a dot plot summarizing the enriched GO terms.
dotplot(ego, showCategory = 15, title = "GO Enrichment Analysis (BP) on Mapped Genes (27 Genes)")
p <- dotplot(ego, showCategory = 15, title = "GO Enrichment Analysis(Top 27 Genes)")
p + theme(axis.text.y = element_text(size = 8))
###########
# Create a bar plot of the GO enrichment results
p <- barplot(ego, showCategory = 15, title = "GO Enrichment Analysis(Top 27 Genes)")

# Adjust the text size if needed
p + theme(axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8))
############
# Load required package (if not already loaded)
library(enrichplot)

# Create a gene network plot (concept network plot) from your GO enrichment results stored in 'ego'.
# 'showCategory' controls how many terms are displayed — here, we show 5 for clarity.
cnet_plot <- cnetplot(ego, 
                      showCategory = 5,
                      circular = FALSE,    # set to TRUE for a circular layout if preferred
                      colorEdge = TRUE) + 
  ggtitle("   Gene Network Plot for GO Enrichment (BP)")

# Print the plot
print(cnet_plot)
#################
library(clusterProfiler)
library(enrichplot)

# Compute the pairwise similarity for the enriched GO terms.
ego2 <- pairwise_termsim(ego)

# Now create the enrichment map plot using the updated object.
emap <- emapplot(ego2, showCategory = 15, layout = "kk") +
  ggtitle("Enrichment Map of GO Terms (BP)") +
  theme_minimal()

print(emap)
#####################
library(clusterProfiler)
library(enrichplot)

# Compute the pairwise similarity for the enriched GO terms.
ego2 <- pairwise_termsim(ego)

# Create the enrichment map plot and remove x and y axis labels and points
emap <- emapplot(ego2, showCategory = 15, layout = "kk") +
  ggtitle("Enrichment Map of GO Terms (BP)") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),   # Remove axis labels
    axis.text  = element_blank(),   # Remove the axis text
    axis.ticks = element_blank(),   # Remove axis ticks
    panel.grid = element_blank()    # Remove grid lines if desired
  )

print(emap)

#######################
# Example SHAP values for one observation (replace these with your actual values)
shap_values <- c(
  "Feature_A" = 0.25,    # Positive contribution (pushes prediction higher)
  "Feature_B" = -0.30,   # Negative contribution (pushes prediction lower)
  "Feature_C" = 0.10,    # Positive contribution
  "Feature_D" = -0.15    # Negative contribution
)

# Base (expected mean) value from your model
base_value <- 0.0412

# Calculate the final prediction as the sum of the base value and the SHAP contributions
predicted_value <- base_value + sum(shap_values)

# Display the base value, sum of contributions, and final prediction
cat("Base value:", base_value, "\n")
cat("Sum of SHAP values:", sum(shap_values), "\n")
cat("Final prediction:", predicted_value, "\n\n")

# Check which features contribute positively and negatively
positive_contributors <- shap_values[shap_values > 0]
negative_contributors <- shap_values[shap_values < 0]

cat("Features with positive contributions (push prediction higher):\n")
print(positive_contributors)

cat("\nFeatures with negative contributions (push prediction lower):\n")
print(negative_contributors)

# Compare the final prediction to the base value
# According to our framework:
#  - If predicted_value > base_value then the positive shifts lead to evasion.
#  - If predicted_value < base_value then the negative shifts lead to non-evasion.

if(predicted_value > base_value) {
  cat("\nResult: The final prediction is higher than the base value. \n")
  cat("This means that the positive contributions (which push the prediction up) are driving the model toward Evasion.\n")
} else if(predicted_value < base_value) {
  cat("\nResult: The final prediction is lower than the base value. \n")
  cat("This indicates that the negative contributions (which push the prediction down) are driving the model toward Non-Evasion.\n")
} else {
  cat("\nResult: The final prediction is equal to the base value, so there is no directional shift.\n")
}
