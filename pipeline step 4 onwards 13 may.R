# =============================================================================
# 0. Setup: load libraries
# =============================================================================
library(mixOmics)
library(caret)
library(BiocParallel)

# =============================================================================
# 1. Load preprocessed data (X blocks, Y, design if you saved it)
# =============================================================================
load("diablo_pipeline_data.RData")    # loads X, Y, diablo_model (old), design (optional)

# =============================================================================
# 2. Sanity checks on X and Y
# =============================================================================
print(lapply(X, dim))   # should be n_samples × n_features in each block
print(length(Y))        # should equal n_samples
stopifnot(
  identical(rownames(X$proteomics),    names(Y)),
  identical(rownames(X$methylation),   names(Y)),
  identical(rownames(X$transcriptomics), names(Y))
)
####
# Which samples are in X$proteomics but not in Y?
setdiff(rownames(X$proteomics), names(Y))
# Which samples are in Y but not in X$proteomics?
setdiff(names(Y), rownames(X$proteomics))

# =============================================================================
# 3. Remove near-zero-variance features in transcriptomics
# =============================================================================
nzv <- nearZeroVar(X$transcriptomics)
if(length(nzv) > 0) {
  X$transcriptomics <- X$transcriptomics[, -nzv, drop = FALSE]
}

# =============================================================================
# 4. Re-build design matrix
# =============================================================================
design <- matrix(0.1,
                 nrow = length(X),
                 ncol = length(X),
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

# =============================================================================
# 4a. Train/Test Split (hold-out 20% for final evaluation)
# =============================================================================
set.seed(42)
train_idx <- createDataPartition(Y, p = 0.8, list = FALSE)
X_train <- lapply(X, `[`, train_idx, , drop = FALSE)
Y_train <- Y[train_idx]
X_test  <- lapply(X, `[-`, train_idx, , drop = FALSE)
Y_test  <- Y[-train_idx]

# Optional: re-check dimensions on training set
print(lapply(X_train, dim))
print(length(Y_train))

# =============================================================================
# 5. Hyperparameter tuning on training set
# =============================================================================
bp_param <- SnowParam(workers = 3, type = "SOCK")

tune.res <- tune.block.splsda(
  X           = X_train,
  Y           = Y_train,
  ncomp       = 2,
  test.keepX  = list(
    proteomics     = seq(5, 50, by = 5),
    methylation    = seq(5, 50, by = 5),
    transcriptomics = seq(5, 50, by = 5)
  ),
  design      = design,
  validation  = "Mfold",
  folds       = 5,
  nrepeat     = 5,
  BPPARAM     = bp_param
)

list.keepX <- tune.res$choice.keepX
print(list.keepX)

# =============================================================================
# 6. Build tuned DIABLO model on training set
# =============================================================================
diablo_model_tuned <- block.splsda(
  X      = X_train,
  Y      = Y_train,
  ncomp  = 2,
  design = design,
  keepX  = list.keepX
)

# =============================================================================
# 7. Performance assessment on training set
# =============================================================================
perf.res <- perf(
  diablo_model_tuned,
  validation  = "Mfold",
  folds       = 5,
  nrepeat     = 10,
  progressBar = TRUE,
  BPPARAM     = bp_param
)
plot(perf.res, main = "Tuned DIABLO Performance (Train Set)")

# =============================================================================
# 8. Extract selected features using selectVar() and loadings
# =============================================================================
extract_features_from_loadings <- function(mod, block, comp) {
  L <- mod$loadings[[block]]
  rownames(L)[L[, comp] != 0]
}

extract_feats_selectVar <- function(mod, block, comp) {
  selectVar(mod, block = block, comp = comp)$value$names
}

blocks <- names(X_train)
# initialize lists
prot_feats_var <- prot_feats_loadings <- NULL
meth_feats_var <- meth_feats_loadings <- NULL
rna_feats_var  <- rna_feats_loadings  <- NULL

# loop blocks for brevity
for(bk in blocks) {
  var_feats <- unique(c(
    extract_feats_selectVar(diablo_model_tuned, bk, 1),
    extract_feats_selectVar(diablo_model_tuned, bk, 2)
  ))
  load_feats <- unique(c(
    extract_features_from_loadings(diablo_model_tuned, bk, 1),
    extract_features_from_loadings(diablo_model_tuned, bk, 2)
  ))
  assign(paste0(bk, "_feats_var"), var_feats)
  assign(paste0(bk, "_feats_loadings"), load_feats)
  cat(sprintf("%s (selectVar): %d features\n", bk, length(var_feats)))
  cat(sprintf("%s (loadings): %d features\n", bk, length(load_feats)))
}

# =============================================================================
# 9. Component-Specific Feature Extraction
# =============================================================================
extract_selected_by_component <- function(model) {
  blocks <- names(model$loadings)
  ncomp  <- model$ncomp[1]
  res    <- list()
  for(bk in blocks) {
    L <- model$loadings[[bk]]
    comps <- lapply(1:ncomp, function(comp) rownames(L)[L[, comp] != 0])
    names(comps) <- paste0("Comp", 1:ncomp)
    res[[bk]] <- comps
  }
  res
}
features_by_comp <- extract_selected_by_component(diablo_model_tuned)
for(bk in names(features_by_comp)) {
  cat(sprintf("\nBlock: %s\n", bk))
  for(comp in names(features_by_comp[[bk]])) {
    feats <- features_by_comp[[bk]][[comp]]
    cat(sprintf("  %s: %d features\n", comp, length(feats)))
  }
}

# =============================================================================
# 10. Final evaluation on held-out test set
# =============================================================================
preds <- predict(diablo_model_tuned, X_test)$class$max.dist
cat("Test-set confusion matrix:\n")
print(confusionMatrix(preds, Y_test))
