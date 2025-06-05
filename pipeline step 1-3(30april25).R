# =============================================================================
# Step 1: Data Loading and Preliminary Exploration
# =============================================================================

# Load necessary libraries for data handling and exploration
library(dplyr)    # For data manipulation
library(readr)    # For efficient CSV reading
library(ggplot2)  # For quick plotting if needed

# --- Specify File Paths ---
# Use forward slashes (or double backslashes) in the file paths
proteomics_file      <- "C:/Users/KIIT/Downloads/multiomics analysis/multiomics/aggregated_proteomics.csv"
methylation_file     <- "C:/Users/KIIT/Downloads/multiomics analysis/multiomics/gene_level_methylation.csv"
transcriptomics_file <- "C:/Users/KIIT/Downloads/multiomics analysis/multiomics/Tumornames.csv"

# --- Load the Data ---
# Read the files as tibbles
proteomics_data <- read_csv(proteomics_file)
methylation_data <- read_csv(methylation_file)
transcriptomics_data <- read_csv(transcriptomics_file)

# --- Preprocess to Set Gene Names as Row Names ---

# For proteomics and transcriptomics, gene names are in the first column.
# Convert to a data.frame and then set the rownames.
proteomics_data <- as.data.frame(proteomics_data)
rownames(proteomics_data) <- proteomics_data[, 1]
proteomics_data <- proteomics_data[, -1]

transcriptomics_data <- as.data.frame(transcriptomics_data)
rownames(transcriptomics_data) <- transcriptomics_data[, 1]
transcriptomics_data <- transcriptomics_data[, -1]

# For methylation, gene names are in the second column.
# Remove the first column, then check for missing gene names in what becomes the first column.
methylation_data <- as.data.frame(methylation_data)
methylation_data <- methylation_data[, -1]  # Remove the extraneous first column

# Check for missing gene names and remove rows with missing values
missing_gene_flag <- is.na(methylation_data[, 1])
if (any(missing_gene_flag)) {
  cat("Found", sum(missing_gene_flag), "rows with missing gene names in methylation data. Removing these rows...\n")
  methylation_data <- methylation_data[!missing_gene_flag, ]
}

# Set the gene names as rownames and remove that column
rownames(methylation_data) <- methylation_data[, 1]
methylation_data <- methylation_data[, -1]
##########################
head(rownames(methylation_data))
# and
dimnames(methylation_data)

##############################
# --- Quick Overview of the Data ---

# Print the dimensions of each dataset
cat("Proteomics data dimensions: ", dim(proteomics_data), "\n")
cat("Methylation data dimensions: ", dim(methylation_data), "\n")
cat("Transcriptomics data dimensions: ", dim(transcriptomics_data), "\n")

# Display the first few rows for a visual check
cat("\n--- Proteomics Data Preview ---\n")
print(head(proteomics_data))
cat("\n--- Methylation Data Preview ---\n")
print(head(methylation_data))
cat("\n--- Transcriptomics Data Preview ---\n")
print(head(transcriptomics_data))

# --- Check Column Names and Structure ---
# This is helpful to see how sample IDs are stored.
cat("\nProteomics Data Column Names:\n")
print(colnames(proteomics_data))
cat("\nMethylation Data Column Names:\n")
print(colnames(methylation_data))
cat("\nTranscriptomics Data Column Names:\n")
print(colnames(transcriptomics_data))
######################
# --- Standardize Sample IDs Across Datasets ---

# Function to standardize sample names:
standardize_sample_id <- function(x) {
  # If the IDs contain dots, replace them with hyphens:
  s <- gsub("\\.", "-", x)
  # Split by hyphen and extract the first three segments:
  segments <- unlist(strsplit(s, "-"))
  if (length(segments) >= 3) {
    standardized <- paste(segments[1:3], collapse = "-")
  } else {
    standardized <- s  # fallback if not enough segments
  }
  return(standardized)
}

# Apply standardization to proteomics sample names:
proteomics_std <- sapply(colnames(proteomics_data), standardize_sample_id)
colnames(proteomics_data) <- proteomics_std

# For methylation and transcriptomics, the IDs appear to be already hyphen-separated.
# We'll extract the first three segments for consistency.
extract_three_segments <- function(x) {
  segments <- unlist(strsplit(x, "-"))
  if (length(segments) >= 3) {
    return(paste(segments[1:3], collapse = "-"))
  } else {
    return(x)
  }
}

methylation_std <- sapply(colnames(methylation_data), extract_three_segments)
colnames(methylation_data) <- methylation_std

transcriptomics_std <- sapply(colnames(transcriptomics_data), extract_three_segments)
colnames(transcriptomics_data) <- transcriptomics_std

# Check standardized sample IDs
cat("\nStandardized Proteomics Sample IDs (first 10):\n")
print(head(colnames(proteomics_data), 10))
cat("\nStandardized Methylation Sample IDs (first 10):\n")
print(head(colnames(methylation_data), 10))
cat("\nStandardized Transcriptomics Sample IDs (first 10):\n")
print(head(colnames(transcriptomics_data), 10))

# Determine common samples across all datasets:
common_samples <- Reduce(intersect, list(colnames(proteomics_data), 
                                         colnames(methylation_data), 
                                         colnames(transcriptomics_data)))
cat("\nNumber of common samples across datasets: ", length(common_samples), "\n")
#################
# --- Align Data Matrices Based on Common Samples ---

# Subset and reorder each dataset using the common sample IDs
proteomics_aligned   <- proteomics_data[, common_samples, drop = FALSE]
methylation_aligned  <- methylation_data[, common_samples, drop = FALSE]
transcriptomics_aligned <- transcriptomics_data[, common_samples, drop = FALSE]

# Check dimensions after alignment:
cat("Aligned Proteomics dimensions: ", dim(proteomics_aligned), "\n")
cat("Aligned Methylation dimensions: ", dim(methylation_aligned), "\n")
cat("Aligned Transcriptomics dimensions: ", dim(transcriptomics_aligned), "\n")

###########################

# --- Summary Statistics ---
cat("\n--- Summary Statistics ---\n")
cat("\nProteomics Data Summary:\n")
print(summary(proteomics_data))
cat("\nMethylation Data Summary:\n")
print(summary(methylation_data))
cat("\nTranscriptomics Data Summary:\n")
print(summary(transcriptomics_data))
###########
# Preview the first 10 rows to inspect the data
cat("\n--- Transcriptomics Data Preview (first 10 rows) ---\n")
print(head(transcriptomics_data, 10))
# Check the distribution for a specific gene (replace 'TSPAN6' with another gene if desired)
cat("\n--- Expression Values for TSPAN6 ---\n")
print(transcriptomics_data["TSPAN6", ])

# Plot density for one sample
library(ggplot2)
df <- data.frame(Expression = transcriptomics_data[, "TCGA-CV-7423"])
ggplot(df, aes(x = Expression)) + 
  geom_density(fill = "skyblue", alpha = 0.5) +
  ggtitle("Density Plot for Sample TCGA-CV-7423")


######################

# --- Optional: Check for Missing Values ---
cat("\nNumber of missing values in Proteomics data: ", sum(is.na(proteomics_data)), "\n")
cat("Number of missing values in Methylation data: ", sum(is.na(methylation_data)), "\n")
cat("Number of missing values in Transcriptomics data: ", sum(is.na(transcriptomics_data)), "\n")

# =============================================================================
# End of Data Loading and Exploration
# =============================================================================

# Once you inspect the printed outputs and verify that the gene names 
# and sample information are correct, you can proceed with further analyses.
####################
library(ggplot2)

# Plot density for the first sample in Proteomics
sample_prot <- colnames(proteomics_aligned)[1]
df_prot <- data.frame(Value = proteomics_aligned[, sample_prot])
ggplot(df_prot, aes(x = Value)) +
  geom_density(fill = "blue", alpha = 0.5) +
  ggtitle(paste("Proteomics - Density Distribution for", sample_prot))

# Plot density for the first sample in Methylation
sample_meth <- colnames(methylation_aligned)[1]
df_meth <- data.frame(Value = methylation_aligned[, sample_meth])
ggplot(df_meth, aes(x = Value)) +
  geom_density(fill = "green", alpha = 0.5) +
  ggtitle(paste("Methylation - Density Distribution for", sample_meth))

# Plot density for the first sample in Transcriptomics
sample_trans <- colnames(transcriptomics_aligned)[1]
df_trans <- data.frame(Value = transcriptomics_aligned[, sample_trans])
ggplot(df_trans, aes(x = Value)) +
  geom_density(fill = "red", alpha = 0.5) +
  ggtitle(paste("Transcriptomics - Density Distribution for", sample_trans))
#################
# --------------------------
# Imputation Option B
# --------------------------
# Define a function for median imputation per gene
impute_median <- function(mat) {
  t(apply(mat, 1, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    return(x)
  }))
}

# For Proteomics: Replace missing values with the median for each gene
proteomics_imputed_median <- impute_median(proteomics_aligned)

# For Methylation: Replace missing values with the median for each gene
methylation_imputed <- impute_median(methylation_aligned)

# Transcriptomics remains unchanged
transcriptomics_imputed <- transcriptomics_aligned

# Check missing values post imputation:
cat("Missing values in Proteomics (median imputation) after imputation: ", sum(is.na(proteomics_imputed_median)), "\n")
cat("Missing values in Methylation after imputation: ", sum(is.na(methylation_imputed)), "\n")
cat("Missing values in Transcriptomics: ", sum(is.na(transcriptomics_imputed)), "\n")
#######################
rownames(methylation_aligned)      # should be gene IDs
rownames(methylation_imputed)      # probably NULL or numeric now

##########################

# Check for rows in proteomics_aligned where all values are missing
all_missing <- apply(proteomics_aligned, 1, function(x) all(is.na(x)))
num_all_missing <- sum(all_missing)
cat("Number of rows with all missing values in Proteomics data:", num_all_missing, "\n")

##
# Remove rows where all values are missing
proteomics_filtered <- proteomics_aligned[!all_missing, ]
cat("Dimensions after removing rows with all missing values:", dim(proteomics_filtered), "\n")
###
# Define a function for row-wise median imputation
impute_median <- function(mat) {
  t(apply(mat, 1, function(x) {
    # Calculate the median of the non-missing values
    med <- median(x, na.rm = TRUE)
    # Replace the remaining NA's with the median
    x[is.na(x)] <- med
    return(x)
  }))
}

# Apply median imputation on the filtered proteomics data
proteomics_imputed <- impute_median(proteomics_filtered)

# Check that there are no missing values left
cat("Missing values in Proteomics after median imputation: ", sum(is.na(proteomics_imputed)), "\n")
#######################################
# 1. Summary statistics for the imputed proteomics data
cat("Summary of proteomics_imputed data:\n")
print(summary(proteomics_imputed))

# 2. Compare density distributions for a few proteins
# Let's inspect the first three proteins for this example:
proteins_to_plot <- rownames(proteomics_imputed)[1:3]

# Load ggplot2 for prettier plots (if not already loaded)
library(ggplot2)

# For each protein, plot density curves for:
# - "Original (non-NA)" values from proteomics_filtered (the data before imputation)
# - "Imputed" values from proteomics_imputed (which include the median imputed values)
for (prot in proteins_to_plot) {
  # Extract the values for this protein from the filtered (pre-imputed) data
  original_values <- as.numeric(proteomics_filtered[prot, ])
  # Remove NAs (only the observed values)
  original_values <- original_values[!is.na(original_values)]
  
  # Extract all values (observed plus imputed) from the imputed data
  imputed_values <- as.numeric(proteomics_imputed[prot, ])
  
  # Combine these into a single data frame for plotting
  df <- data.frame(
    Value = c(original_values, imputed_values),
    Source = factor(c(rep("Original (non-NA)", length(original_values)),
                      rep("Imputed", length(imputed_values))))
  )
  
  p <- ggplot(df, aes(x = Value, fill = Source)) +
    geom_density(alpha = 0.5) +
    ggtitle(paste("Density Comparison for", prot)) +
    theme_minimal()
  
  print(p)
}

#########################################
# Load ggplot2 for plotting (if not already loaded)
library(ggplot2)

# Select a few genes (features) to compare - here, we use the first three rows for illustration
genes_to_plot <- rownames(methylation_imputed)[1:3]

# Loop over the selected genes and create a density plot for each
for (gene in genes_to_plot) {
  # Extract original values from the un-imputed methylation data for the gene
  original_values <- as.numeric(methylation_aligned[gene, ])
  # Remove NAs to get only observed values
  original_values <- original_values[!is.na(original_values)]
  
  # Extract the complete set of values (observed plus imputed) for the gene
  imputed_values <- as.numeric(methylation_imputed[gene, ])
  
  # Create a combined data frame for plotting the two groups
  df <- data.frame(
    Value = c(original_values, imputed_values),
    Source = factor(c(rep("Original (non-NA)", length(original_values)),
                      rep("Imputed", length(imputed_values))))
  )
  
  # Create the density plot
  p <- ggplot(df, aes(x = Value, fill = Source)) +
    geom_density(alpha = 0.5) +
    ggtitle(paste("Density Comparison for", gene)) +
    theme_minimal() +
    xlab("Methylation Beta Value") +
    ylab("Density")
  
  print(p)
}
###############################
# For proteomics data
row_means_prot <- apply(proteomics_imputed, 1, mean, na.rm = TRUE)
row_sds_prot   <- apply(proteomics_imputed, 1, sd, na.rm = TRUE)

cat("Proteomics row means summary:\n")
print(summary(row_means_prot))
cat("\nProteomics row standard deviations summary:\n")
print(summary(row_sds_prot))
##
# For methylation data
row_means_meth <- apply(methylation_imputed, 1, mean, na.rm = TRUE)
row_sds_meth   <- apply(methylation_imputed, 1, sd, na.rm = TRUE)

cat("Methylation row means summary:\n")
print(summary(row_means_meth))
cat("\nMethylation row standard deviations summary:\n")
print(summary(row_sds_meth))


# For transcriptomics data
row_means_rna <- apply(transcriptomics_imputed, 1, mean, na.rm = TRUE)
row_sds_rna   <- apply(transcriptomics_imputed, 1, sd, na.rm = TRUE)

cat("Transcriptomics row means summary:\n")
print(summary(row_means_rna))
cat("\nTranscriptomics row standard deviations summary:\n")
print(summary(row_sds_rna))
##############################
# Define a function for row-wise Z‑score scaling (each gene across samples)
scale_data <- function(mat) {
  t(scale(t(mat)))
}

# Apply scaling to each omics matrix
proteomics_scaled    <- scale_data(proteomics_imputed)
methylation_scaled   <- scale_data(methylation_imputed)
transcriptomics_scaled <- scale_data(transcriptomics_imputed)

# Verify the scaling by re-checking row means and standard deviations
prot_means <- apply(proteomics_scaled, 1, mean, na.rm = TRUE)
prot_sds   <- apply(proteomics_scaled, 1, sd, na.rm = TRUE)
cat("Proteomics scaled: Row means summary:\n")
print(summary(prot_means))
cat("Proteomics scaled: Row standard deviations summary:\n")
print(summary(prot_sds))

meth_means <- apply(methylation_scaled, 1, mean, na.rm = TRUE)
meth_sds   <- apply(methylation_scaled, 1, sd, na.rm = TRUE)
cat("Methylation scaled: Row means summary:\n")
print(summary(meth_means))
cat("Methylation scaled: Row standard deviations summary:\n")
print(summary(meth_sds))

rna_means <- apply(transcriptomics_scaled, 1, mean, na.rm = TRUE)
rna_sds   <- apply(transcriptomics_scaled, 1, sd, na.rm = TRUE)
cat("Transcriptomics scaled: Row means summary:\n")
print(summary(rna_means))
cat("Transcriptomics scaled: Row standard deviations summary:\n")
print(summary(rna_sds))
##########################################
rownames(methylation_scaled)

####################################################
# Identify rows with non-zero variance
nonzero_var_idx <- apply(transcriptomics_imputed, 1, function(x) {
  var(x, na.rm = TRUE) > 0
})

# Filter out rows that have zero variance
transcriptomics_nonzero <- transcriptomics_imputed[nonzero_var_idx, ]
cat("Dimensions before filtering zero-variance genes: ", dim(transcriptomics_imputed), "\n")
cat("Dimensions after filtering zero-variance genes: ", dim(transcriptomics_nonzero), "\n")

# Now apply the Z-score scaling function on the filtered data
scale_data <- function(mat) {
  t(scale(t(mat)))
}

transcriptomics_scaled <- scale_data(transcriptomics_nonzero)

# Verify the scaling after filtering
rna_means <- apply(transcriptomics_scaled, 1, mean, na.rm = TRUE)
rna_sds   <- apply(transcriptomics_scaled, 1, sd, na.rm = TRUE)
cat("Transcriptomics scaled: Row means summary:\n")
print(summary(rna_means))
cat("Transcriptomics scaled: Row standard deviations summary:\n")
print(summary(rna_sds))
#######################################
library(Seurat)

# Create a Seurat object from your transcriptomics data.
# 'transcriptomics_scaled' is your pre-processed, Z-score normalized data matrix.
seurat_object <- CreateSeuratObject(counts = transcriptomics_scaled)

# Define your curated list of immune evasion genes.
immune_evasion_genes <- c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "CD244", 
                          "CD160", "IRF4", "BATF", "NFATC1", "FOXP3", "IL2RA", 
                          "TNFRSF18", "IKZF4", "CD80", "CD86", "CD274", 
                          "PDCD1LG2", "LGALS9")

# Ensure the default assay is set to "RNA"
DefaultAssay(seurat_object) <- "RNA"

# Assign your scaled transcriptomics data to the 'data' slot using SetAssayData.
seurat_object <- SetAssayData(
  object = seurat_object,
  assay = "RNA",
  slot = "data",
  new.data = as.matrix(transcriptomics_scaled)
)

# Calculate the immune evasion module score using AddModuleScore.
# This function will compute a score based on the mean expression
# of the genes in the provided list.
seurat_object <- AddModuleScore(
  object = seurat_object,
  features = list(ImmuneEvasion = immune_evasion_genes),
  name = "ImmuneEvasionScore"
)

# Extract the immune evasion scores from the Seurat metadata.
immune_scores <- seurat_object@meta.data$ImmuneEvasionScore1

# Plot the distribution of the immune evasion scores.
hist(immune_scores, breaks = 30, main = "Immune Evasion Score Distribution", xlab = "Score")
####################
library(ggplot2)

# Extract meta data from the Seurat object
meta_df <- seurat_object@meta.data

# Create a violin plot for the immune evasion scores across all samples
ggplot(meta_df, aes(x = "", y = ImmuneEvasionScore1)) +
  geom_violin(fill = "skyblue", alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(x = NULL, y = "Immune Evasion Score",
       title = "Violin Plot of Immune Evasion Scores") +
  theme_minimal()
###############################

#########################
# Binarize the scores based on the median value.
threshold <- median(immune_scores, na.rm = TRUE)
Y <- ifelse(immune_scores >= threshold, "Evasive", "Non-Evasive")
Y <- factor(Y)
print(table(Y))

# Now Y can be used as your phenotype for DIABLO integration.
###############################
library(ggplot2)

# Add your binarized scores as a new metadata column in the Seurat object.
seurat_object@meta.data$ImmuneGroup <- Y

# Extract the metadata into a data frame for plotting.
meta_df <- seurat_object@meta.data

# Create the grouped violin plot.
ggplot(meta_df, aes(x = ImmuneGroup, y = ImmuneEvasionScore1, fill = ImmuneGroup)) +
  geom_violin(trim = FALSE, adjust = 1.0, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  labs(x = "Immune Evasion Status", 
       y = "Immune Evasion Score", 
       title = "Violin Plot of Immune Evasion Scores by Group") +
  theme_minimal() +
  scale_fill_manual(values = c("Non-Evasive" = "skyblue", "Evasive" = "salmon"))
#######################
#–– Check methylation before X creation
stopifnot(!is.null(rownames(methylation_scaled)))
stopifnot(is.character(rownames(methylation_scaled)))
#####################################
library(mixOmics)

# Step 1: Assemble the omics blocks into a list.
# Transpose each scaled matrix so that samples (rows) and features (columns) are as expected.
X <- list(
  proteomics = t(proteomics_scaled),
  methylation = t(methylation_scaled),
  transcriptomics = t(transcriptomics_scaled)
)
####################
head(colnames(X$methylation)) 
# 1a. Dimensions
lapply(X, dim)
# should all report:  n_samples × n_features

length(Y)
# should equal n_samples

# 1b. Names
# If your X blocks have rownames (samples), make sure they match Y
all(rownames(X$proteomics)   == names(Y))
all(rownames(X$methylation)  == names(Y))
all(rownames(X$transcriptomics) == names(Y))

########################
# Optional: Verify the dimensions of each block. They should all have 335 rows (samples).
lapply(X, dim)

# Step 2: Outcome variable.
# Y has been created as a factor (Evasive vs. Non-Evasive) for 335 samples.
print(table(Y))

# Step 3: Define the design matrix.
# Off-diagonals set to 0.1 and zeros on the diagonal.
design <- matrix(0.1, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

# Optional: Inspect the design matrix.
print(design)

# Step 4: Fit the DIABLO model.
# Fit with 2 components based on your design.
diablo_model <- block.splsda(X, Y, ncomp = 2, design = design)

# Quick summary of the model:
print(diablo_model)

# Visualize the sample projections on the latent components:
plotIndiv(diablo_model, legend = TRUE, title = "DIABLO Sample Projections")
#######
plotIndiv(diablo_model, legend = TRUE, title = "DIABLO Sample Projections",
          ind.names = FALSE, ellipse = TRUE)
############################
plotIndiv(diablo_model, 
          legend = TRUE, 
          title = "DIABLO Sample Projections",
          ind.names = FALSE,      # Do not print sample IDs.
          ellipse = TRUE,         # Add 95% confidence ellipses.
          style = "ggplot2",      # Use ggplot2 style for a cleaner look.
          pch = c(21, 24),        # Use different point shapes (e.g., circle & triangle).
          col = c("salmon", "skyblue"),  # Set contrasting colors for groups.
          cex = 3)                # Increase point size for better visibility.

#############
# Optional: Explore the variable importance in one block.
plotLoadings(diablo_model, comp = 1, block = "proteomics", title = "Proteomics Loadings on Component 1")
##############
# Save key objects to a file
save(X, Y, diablo_model, file = "diablo_pipeline_data.RData")

# Optionally, if you want to save additional objects, e.g., your design matrix:
# save(X, Y, diablo_model, design, file = "diablo_pipeline_data.RData")
load("diablo_pipeline_data.RData")
#############
# Load the saved data
load("diablo_pipeline_data.RData")  # This loads X, Y, diablo_model (and optionally design)

colnames(X$methylation)[1:10]


# Check dimensions of X (which is usually a list of datasets in DIABLO)
lapply(X, dim)

# Check dimensions or structure of Y
dim(Y)  # if Y is a matrix or data frame
# or
length(Y)  # if Y is a factor or vector

# Check the structure of the diablo_model (to confirm it's loaded properly)
str(diablo_model)

###################stop here ,check pipline 1.5
library(caret)

# Identify near zero variance features in the transcriptomics block
nzv_idx <- nearZeroVar(X$transcriptomics)

# Print how many near-zero variance features are detected
cat("Number of near zero variance features in transcriptomics:", length(nzv_idx), "\n")
# Check the total number of features (columns) in the transcriptomics data
total_features <- ncol(X$transcriptomics)
cat("Number of features in transcriptomics data:", total_features, "\n")

# Remove these features if any are found
if(length(nzv_idx) > 0) {
  cat("Removing near zero variance features in transcriptomics with indices:", nzv_idx, "\n")
  X$transcriptomics <- X$transcriptomics[, -nzv_idx, drop = FALSE]
} else {
  cat("No near zero variance features detected in transcriptomics.\n")
}

# Check and print how many features remain in transcriptomics
cat("Transcriptomics block now has", ncol(X$transcriptomics), "features.\n")
########################
# Re-create the design matrix (if needed)
design <- matrix(0.1, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

# Re-fit the DIABLO model on the updated data
diablo_model_filtered <- block.splsda(X = X, Y = Y, ncomp = 2, design = design)
#####
perf_diablo <- perf(diablo_model_filtered,
                    validation = "Mfold",
                    folds = 5,
                    nrepeat = 10,
                    near.zero.var = TRUE,
                    progressBar = TRUE)

print(perf_diablo)
names(perf_diablo$error.rate)

plot(perf_diablo$error.rate$proteomics,
     main = "DIABLO CV Error Rates: Proteomics",
     xlab = "Component",
     ylab = "Error Rate")
###############refitted diablo due to some existing variance 
par(mfrow = c(1, 3))
plot(perf_diablo$error.rate$proteomics,
     main = "Proteomics",
     xlab = "Component",
     ylab = "Error Rate")
plot(perf_diablo$error.rate$methylation,
     main = "Methylation",
     xlab = "Component",
     ylab = "Error Rate")
plot(perf_diablo$error.rate$transcriptomics,
     main = "Transcriptomics",
     xlab = "Component",
     ylab = "Error Rate")
par(mfrow = c(1, 1))  # Reset to one plot per window if needed.
#############################
#saving for safety
saveRDS(X, file = "X_preTuning.rds")
saveRDS(Y, file = "Y_preTuning.rds")
saveRDS(design, file = "design_preTuning.rds")
######################
# Load required libraries
library(mixOmics)
# Load data
X <- readRDS("X_preTuning.rds")
Y <- readRDS("Y_preTuning.rds")
design <- readRDS("design_preTuning.rds")
#############
lapply(X, dim)
length(Y)
table(Y)        # to see class distribution
is.factor(Y)    # to confirm if it's a factor
dim(design)
design

##################11
# Load BiocParallel
library(BiocParallel)

# Set up SnowParam for parallel processing (Windows)
bp_param <- SnowParam(workers = 3, type = "SOCK")

# 4.1 Hyperparameter Tuning with tune.block.splsda()
tune.res <- tune.block.splsda(
  X = X, 
  Y = Y, 
  ncomp = 2,
  test.keepX = list(
    proteomics = seq(5, 50, by = 5),
    methylation = seq(5, 50, by = 5),
    transcriptomics = seq(5, 50, by = 5)
  ),
  design = design, 
  validation = "Mfold", 
  folds = 5, 
  nrepeat = 10,
  BPPARAM = bp_param
)

# Extract the optimal number of features per block
list.keepX <- tune.res$choice.keepX
cat("Optimal feature numbers selected per block:\n")
print(list.keepX)

# 4.2 Build a tuned DIABLO model using the optimal keepX values
diablo_model_tuned <- block.splsda(
  X = X, 
  Y = Y, 
  ncomp = 2, 
  design = design,
  keepX = list.keepX
)

# 4.2.1 Assess performance using cross-validation
perf.res <- perf(
  diablo_model_tuned, 
  validation = "Mfold", 
  folds = 5, 
  nrepeat = 10, 
  progressBar = TRUE,
  BPPARAM = bp_param
)

# Plot overall performance results
plot(perf.res, main = "DIABLO Performance Assessment")

# 4.3 Feature Extraction: concatenate features from components 1 & 2

# Proteomics
prot_comp1 <- selectVar(diablo_model_tuned, block = "proteomics", comp = 1)$value$names
prot_comp2 <- selectVar(diablo_model_tuned, block = "proteomics", comp = 2)$value$names
selected_prot_all <- unique(c(prot_comp1, prot_comp2))

# Methylation
meth_comp1 <- selectVar(diablo_model_tuned, block = "methylation", comp = 1)$value$names
meth_comp2 <- selectVar(diablo_model_tuned, block = "methylation", comp = 2)$value$names
selected_meth_all <- unique(c(meth_comp1, meth_comp2))

# Transcriptomics
rna_comp1 <- selectVar(diablo_model_tuned, block = "transcriptomics", comp = 1)$value$names
rna_comp2 <- selectVar(diablo_model_tuned, block = "transcriptomics", comp = 2)$value$names
selected_rna_all <- unique(c(rna_comp1, rna_comp2))

# Display counts
cat("Total features selected in proteomics:\n")
print(length(selected_prot_all))

cat("Total features selected in methylation:\n")
print(length(selected_meth_all))

cat("Total features selected in transcriptomics:\n")
print(length(selected_rna_all))
#######################
# Blocks to extract
blocks <- c("proteomics", "methylation", "transcriptomics")

# Pull selected features from the loadings
selected_features <- lapply(blocks, function(bk) {
  L <- diablo_model_tuned$loadings[[bk]]
  rownames(L)[ apply(L != 0, 1, any) ]
})
names(selected_features) <- blocks

# Report
for(bk in blocks) {
  cat(sprintf("Block %-14s : %4d selected features\n",
              bk, length(selected_features[[bk]])))
  cat(" (first 10):", head(selected_features[[bk]], 10), "\n\n")
}

#############################
# Check the content of keepX and block order
str(list.keepX)
names(X)
# Check column names exist
colnames(X$proteomics)
colnames(X$methylation)
colnames(X$transcriptomics)
###################
# Do the names match, in the same order?
print(names(X))
print(names(list.keepX))
identical(names(X), names(list.keepX))
