# =============================================================================
# DEMO: Fast end-to-end DIABLO test on small feature subsets
# =============================================================================
library(mixOmics)
library(caret)
library(BiocParallel)

# 1. Load data
load("diablo_pipeline_data.RData")   # expects X (list), Y (factor), design (if saved)

# 2. Subset each block for speed
set.seed(42)
demo_X <- list(
  proteomics      = X$proteomics[,   sample(ncol(X$proteomics),    50)],
  methylation     = X$methylation[,  sample(ncol(X$methylation), 1000)],
  transcriptomics = X$transcriptomics[, sample(ncol(X$transcriptomics), 1000)]
)
demo_Y <- Y

# 3. Tiny design matrix
demo_design <- matrix(0.1, nrow = 3, ncol = 3,
                      dimnames = list(names(demo_X), names(demo_X)))
diag(demo_design) <- 0

# 4. Remove near-zero‐variance in transcriptomics
nzv <- nearZeroVar(demo_X$transcriptomics)
if(length(nzv)) demo_X$transcriptomics <- demo_X$transcriptomics[, -nzv]

# 5. Hyperparameter tuning (very narrow grid, light CV)
bp <- SnowParam(workers = 2, type = "SOCK")
demo_tune <- tune.block.splsda(
  X          = demo_X,
  Y          = demo_Y,
  ncomp      = 2,
  test.keepX = list(
    proteomics      = c(5, 10),
    methylation     = c(20, 50),
    transcriptomics = c(20, 50)
  ),
  design     = demo_design,
  validation = "Mfold",
  folds      = 3,
  nrepeat    = 3,
  BPPARAM    = bp
)
cat("Chosen keepX per block:\n"); print(demo_tune$choice.keepX)

# 6. Fit tuned model
demo_mod <- block.splsda(
  X      = demo_X,
  Y      = demo_Y,
  ncomp  = 2,
  design = demo_design,
  keepX  = demo_tune$choice.keepX
)

# 7. Quick performance check
demo_perf <- perf(
  demo_mod,
  validation = "Mfold",
  folds      = 3,
  nrepeat    = 3,
  progressBar= TRUE,
  BPPARAM    = bp
)
plot(demo_perf, main = "Demo DIABLO Performance")

# 8. Extract selected features directly from loadings
blocks <- names(demo_X)
selected <- lapply(blocks, function(bk) {
  L <- demo_mod$loadings[[bk]]
  rownames(L)[ apply(L != 0, 1, any) ]
})
names(selected) <- blocks

# 9. Report
for(bk in blocks) {
  nfeat <- length(selected[[bk]])
  total <- nrow(demo_mod$loadings[[bk]])
  cat(sprintf("Block %-14s: %3d of %4d features selected\n",
              bk, nfeat, total))
  cat("  First few:", head(selected[[bk]], 5), "\n\n")
}
##############
# =============================================================================
# Component-Specific DIABLO Feature Selection Extraction
# =============================================================================

# Define the function for component-specific feature extraction
extract_selected_by_component <- function(model) {
  blocks <- names(model$loadings)
  ncomp <- model$ncomp[1]  # Safely use the first component count
  out_list <- list()
  
  for (bk in blocks) {
    L <- model$loadings[[bk]]
    comp_feats <- list()
    
    # For each component, extract selected features
    for (comp in 1:ncomp) {
      sel <- rownames(L)[L[, comp] != 0]
      comp_feats[[paste0("Comp", comp)]] <- sel
    }
    
    out_list[[bk]] <- comp_feats
  }
  
  return(out_list)
}

# 1. Use your actual trained model for the component-specific extraction
model <- demo_mod  # Replace with your trained DIABLO model

# 2. Perform the component-specific feature extraction
features_by_comp <- extract_selected_by_component(model)

# 3. Report component-specific selected features
blocks <- names(features_by_comp)
for (bk in blocks) {
  cat(sprintf("\nBlock: %s\n", bk))
  for (comp in 1:length(features_by_comp[[bk]])) {
    comp_features <- features_by_comp[[bk]][[paste0("Comp", comp)]]
    cat(sprintf("  Component %d: %d features selected\n", comp, length(comp_features)))
    cat("    First few:", head(comp_features, 5), "\n")
  }
}
