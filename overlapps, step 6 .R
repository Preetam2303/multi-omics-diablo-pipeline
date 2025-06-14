# Load the required data objects
load("lasso_models_results.RData")
load("lasso_selected_features.RData")
# Load enrichment results
load("proteomics_go_enrichment.RData")
load("transcriptomics_enrichment.RData")
load("methylation_kegg_enrichment.RData")
ls()

# Check the names of the blocks
print(names(selected_lasso_features))

# Extract the features for each block
proteomics_features <- selected_lasso_features$proteomics
methylation_features <- selected_lasso_features$methylation
transcriptomics_features <- selected_lasso_features$transcriptomics

# Find common features across all three blocks
common_all_three <- Reduce(intersect, list(proteomics_features, methylation_features, transcriptomics_features))

# Find pairwise overlaps
common_prot_methy <- intersect(proteomics_features, methylation_features)
common_prot_trans <- intersect(proteomics_features, transcriptomics_features)
common_methy_trans <- intersect(methylation_features, transcriptomics_features)

# Print overlap results
cat("Number of features in Proteomics:", length(proteomics_features), "\n")
cat("Number of features in Methylation:", length(methylation_features), "\n")
cat("Number of features in Transcriptomics:", length(transcriptomics_features), "\n")

cat("\nNumber of features common to ALL THREE blocks:", length(common_all_three), "\n")
if(length(common_all_three) > 0){
  print(common_all_three)
} else {
  cat("No features overlap across all three blocks.\n")
}

cat("\nNumber of common features (Proteomics & Methylation):", length(common_prot_methy), "\n")
if(length(common_prot_methy) > 0){
  print(common_prot_methy)
}

cat("\nNumber of common features (Proteomics & Transcriptomics):", length(common_prot_trans), "\n")
if(length(common_prot_trans) > 0){
  print(common_prot_trans)
}

cat("\nNumber of common features (Methylation & Transcriptomics):", length(common_methy_trans), "\n")
if(length(common_methy_trans) > 0){
  print(common_methy_trans)
}

# Optional: visualize overlaps as Venn diagram
library(VennDiagram)
venn.plot <- draw.triple.venn(
  area1 = length(proteomics_features),
  area2 = length(methylation_features),
  area3 = length(transcriptomics_features),
  n12 = length(common_prot_methy),
  n23 = length(common_methy_trans),
  n13 = length(common_prot_trans),
  n123 = length(common_all_three),
  category = c("Proteomics", "Methylation", "Transcriptomics"),
  fill = c("skyblue", "pink1", "lightgreen"),
  lty = "dashed",
  cex = 1.5,
  cat.cex = 1.5,
  euler.d = TRUE,
  scaled = TRUE
)

grid.draw(venn.plot)
#############################################
# Proteomics GO BP pathways
prot_pathways <- prot_go@result$Description

# Transcriptomics GO BP pathways
transcript_go_pathways <- enrich_res@result$Description

# Transcriptomics KEGG pathways
transcript_kegg_pathways <- kegg_res@result$Description

# Methylation KEGG pathways
methy_kegg_pathways <- kegg_methy_res@result$Description
######################
# Overlap at GO BP level (Proteomics & Transcriptomics)
common_go_pathways <- intersect(prot_pathways, transcript_go_pathways)

cat("Number of common GO BP pathways (Proteomics & Transcriptomics):", length(common_go_pathways), "\n")
print(common_go_pathways)

# Overlap at KEGG level (Transcriptomics & Methylation)
common_kegg_pathways <- intersect(transcript_kegg_pathways, methy_kegg_pathways)

cat("Number of common KEGG pathways (Transcriptomics & Methylation):", length(common_kegg_pathways), "\n")
print(common_kegg_pathways)
###################
################################################111111111111111111111111111111
library(clusterProfiler)

# Proteomics genes: convert to Entrez IDs (you already did this)
prot_genes <- selected_lasso_features$proteomics
prot_entrez_df <- convert_to_entrez(prot_genes)

# KEGG enrichment
prot_kegg_res <- enrichKEGG(
  gene         = prot_entrez_df$ENTREZID,
  organism     = 'hsa',
  pvalueCutoff = 0.2
)

# View the results
print(prot_kegg_res)

# Plot top KEGG pathways
if (nrow(prot_kegg_res) > 0) {
  dotplot(prot_kegg_res, showCategory = 10) + ggtitle("Proteomics KEGG Pathway Enrichment")
} else {
  cat("No significant KEGG pathways found for proteomics block.\n")
}

# Save the results
save(prot_kegg_res, file = "proteomics_kegg_enrichment.RData")
############
# Load the KEGG results
load("proteomics_kegg_enrichment.RData")
load("transcriptomics_enrichment.RData")
load("methylation_kegg_enrichment.RData")

# Extract KEGG pathway descriptions
prot_kegg_pathways <- prot_kegg_res@result$Description
transcript_kegg_pathways <- kegg_res@result$Description
methy_kegg_pathways <- kegg_methy_res@result$Description

# Three-way overlap
common_kegg_all <- Reduce(intersect, list(prot_kegg_pathways, transcript_kegg_pathways, methy_kegg_pathways))
cat("Number of common KEGG pathways across all three blocks:", length(common_kegg_all), "\n")
print(common_kegg_all)
#############
library(VennDiagram)

venn.plot_kegg_all <- draw.triple.venn(
  area1 = length(prot_kegg_pathways),
  area2 = length(transcript_kegg_pathways),
  area3 = length(methy_kegg_pathways),
  n12 = length(intersect(prot_kegg_pathways, transcript_kegg_pathways)),
  n23 = length(intersect(transcript_kegg_pathways, methy_kegg_pathways)),
  n13 = length(intersect(prot_kegg_pathways, methy_kegg_pathways)),
  n123 = length(common_kegg_all),
  category = c("Proteomics KEGG", "Transcriptomics KEGG", "Methylation KEGG"),
  fill = c("skyblue", "lightgreen", "pink"),
  lty = "dashed"
)

grid.draw(venn.plot_kegg_all)
######################
library(VennDiagram)

# Create the Venn diagram object
venn.plot_kegg_all <- draw.triple.venn(
  area1 = length(prot_kegg_pathways),
  area2 = length(transcript_kegg_pathways),
  area3 = length(methy_kegg_pathways),
  n12 = length(intersect(prot_kegg_pathways, transcript_kegg_pathways)),
  n23 = length(intersect(transcript_kegg_pathways, methy_kegg_pathways)),
  n13 = length(intersect(prot_kegg_pathways, methy_kegg_pathways)),
  n123 = length(common_kegg_all),
  category = c("Proteomics KEGG", "Transcriptomics KEGG", "Methylation KEGG"),
  fill = c("skyblue", "lightgreen", "pink"),
  lty = "dashed",
  
  # Adjust category label positions and sizes
  cat.cex = 1.2,                # Increase category label size
  cat.fontface = "bold",        # Make labels bold
  cat.pos = c(-20, 20, 180),    # Adjust label positions manually
  cat.dist = c(0.05, 0.05, 0.05),  # Adjust distance from diagram
  cat.just = list(c(0.5, 0), c(0.5, 0), c(0.5, 1)), # Adjust justification
  margin = 0.1                  # Increase margin to prevent cutoff
)

# Draw the modified Venn diagram
grid.draw(venn.plot_kegg_all)
########################################################
# Already done KEGG pathway lists:
# prot_kegg_pathways
# transcript_kegg_pathways
# methy_kegg_pathways

# Compute 2-way overlaps
overlap_prot_trans <- intersect(prot_kegg_pathways, transcript_kegg_pathways)
overlap_prot_methy <- intersect(prot_kegg_pathways, methy_kegg_pathways)
overlap_trans_methy <- intersect(transcript_kegg_pathways, methy_kegg_pathways)

# Show number of overlapping pathways
cat("Proteomics-Transcriptomics overlap:", length(overlap_prot_trans), "\n")
cat("Proteomics-Methylation overlap:", length(overlap_prot_methy), "\n")
cat("Transcriptomics-Methylation overlap:", length(overlap_trans_methy), "\n")
################################
library(ggplot2)

# Create data frame for overlaps
overlap_data <- data.frame(
  Pair = c("Proteomics-Transcriptomics", "Proteomics-Methylation", "Transcriptomics-Methylation"),
  Count = c(length(overlap_prot_trans), length(overlap_prot_methy), length(overlap_trans_methy))
)

# Plot
ggplot(overlap_data, aes(x = Pair, y = Count, fill = Pair)) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_minimal() +
  labs(title = "2-Way KEGG Pathway Overlaps", x = "Pair", y = "Number of Overlapping Pathways") +
  scale_fill_manual(values = c("skyblue", "lightgreen", "pink")) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
#########################
cat("\nProteomics-Transcriptomics KEGG Overlaps:\n")
print(overlap_prot_trans)

cat("\nProteomics-Methylation KEGG Overlaps:\n")
print(overlap_prot_methy)

cat("\nTranscriptomics-Methylation KEGG Overlaps:\n")
print(overlap_trans_methy)
#########################
# 1. Compute pairwise overlaps
prot_transcript_overlap <- intersect(prot_kegg_pathways, transcript_kegg_pathways)
prot_methy_overlap <- intersect(prot_kegg_pathways, methy_kegg_pathways)
transcript_methy_overlap <- intersect(transcript_kegg_pathways, methy_kegg_pathways)

# 2. Compute 3-way overlap
common_kegg_all <- Reduce(intersect, list(prot_kegg_pathways, transcript_kegg_pathways, methy_kegg_pathways))

# 3. Create a summary table
overlap_summary <- data.frame(
  Comparison = c(
    "Proteomics & Transcriptomics",
    "Proteomics & Methylation",
    "Transcriptomics & Methylation",
    "3-way overlap"
  ),
  Overlap_Count = c(
    length(prot_transcript_overlap),
    length(prot_methy_overlap),
    length(transcript_methy_overlap),
    length(common_kegg_all)
  ),
  Pathways = c(
    paste(prot_transcript_overlap, collapse = ", "),
    paste(prot_methy_overlap, collapse = ", "),
    paste(transcript_methy_overlap, collapse = ", "),
    paste(common_kegg_all, collapse = ", ")
  )
)

# 4. Print the table
print(overlap_summary)

# 5. Optionally, save it as CSV for future reports
write.csv(overlap_summary, "KEGG_overlap_summary.csv", row.names = FALSE)
#######################
# List all objects currently in the environment
ls()

length(selected_lasso_features)
head(selected_lasso_features)
#
# Find common features between proteomics and methylation
common_prot_meth <- intersect(selected_lasso_features$proteomics, selected_lasso_features$methylation)
print(common_prot_meth)

# Find common features between proteomics and transcriptomics
common_prot_trans <- intersect(selected_lasso_features$proteomics, selected_lasso_features$transcriptomics)
print(common_prot_trans)

# Find common features between methylation and transcriptomics
common_meth_trans <- intersect(selected_lasso_features$methylation, selected_lasso_features$transcriptomics)
print(common_meth_trans)
