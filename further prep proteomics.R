gene_level_methylation <- readRDS("gene_level_methylation.rds")
output_csv <- "C:/Users/C11_LAB/Downloads/multiomics/gene_level_methylation.csv"
write.csv(gene_level_methylation, file = output_csv, row.names = TRUE)
cat("Data successfully written to", output_csv, "\n")
###################
# Load the CSV file and assign it to 'transcriptomics'
transcriptomics <- read.csv("C:/Users/C11_LAB/Downloads/multiomics/Tumornames.csv", 
                            header = TRUE, 
                            stringsAsFactors = FALSE)

# Inspect the first few rows of the data
head(transcriptomics)
###
gene_level_methylation_csv <- read.csv("C:/Users/C11_LAB/Downloads/multiomics/gene_level_methylation.csv", header = TRUE, stringsAsFactors = FALSE)
head(gene_level_methylation_csv)
###
# Load the proteomics data CSV file into R
proteomics_data <- read.csv("C:/Users/C11_LAB/Downloads/multiomics/proteomics_combined.csv", 
                            header = TRUE, 
                            stringsAsFactors = FALSE, 
                            row.names = 1)

# Inspect the first few rows of the loaded data
head(proteomics_data)
################
# Print all column names for each dataset
cat("Transcriptomics data column names:\n")
print(colnames(transcriptomics))
cat("\nGene-level methylation data column names:\n")
print(colnames(gene_level_methylation_csv))
cat("\nProteomics data column names:\n")
print(colnames(proteomics_data))

#####################################
# 1. Extract sample names from each data set

## Transcriptomics - assuming the first column holds the sample names
transcriptomics_samples <- transcriptomics[[1]]  
# If your transcriptomics file has a different column name,
# for example "SampleID", then you could use:
# transcriptomics_samples <- transcriptomics$SampleID

## Gene-level methylation data -
## Option A: if the first column contains gene identifiers and the remaining columns are samples
gene_level_methylation_samples <- colnames(gene_level_methylation_csv)[-1]
## Option B: if all columns are samples (i.e. there aren’t gene identifiers as a column)
# gene_level_methylation_samples <- colnames(gene_level_methylation_csv)

## Proteomics data - the column names are the sample IDs
proteomics_samples <- colnames(proteomics_data)

# 2. Find the intersection (i.e. common samples) among all three datasets

common_samples <- Reduce(intersect, list(as.character(transcriptomics_samples), 
                                         gene_level_methylation_samples, 
                                         proteomics_samples))
#################

# 3. Check and print the number of common samples
cat("Number of common samples:", length(common_samples), "\n")
print(common_samples)
####################
# --- Standardize Transcriptomics Sample Names ---
# Assume that the first column in 'transcriptomics' contains gene names,
# so the sample names are in all columns except the first.
transcriptomics_samples_full <- colnames(transcriptomics)[-1]
short_transcriptomics <- sapply(transcriptomics_samples_full, function(x) {
  parts <- unlist(strsplit(x, "\\."))
  # Keep only the first 4 fields (e.g., "TCGA.CV.7423.11A" becomes "TCGA.CV.7423.11A")
  paste(parts[1:4], collapse = ".")
})
cat("Short transcriptomics sample names (first few):\n")
print(head(short_transcriptomics))


# --- Standardize Proteomics Sample Names ---
# For proteomics, the column names already represent sample IDs.
proteomics_samples_full <- colnames(proteomics_data)
short_proteomics <- sapply(proteomics_samples_full, function(x) {
  parts <- unlist(strsplit(x, "\\."))
  # Keep only the first 4 fields (e.g., "TCGA.F7.8489.01A.21.A45L.20" becomes "TCGA.F7.8489.01A")
  paste(parts[1:4], collapse = ".")
})
cat("\nShort proteomics sample names (first few):\n")
print(head(short_proteomics))


# --- Extract Gene-level Methylation Sample Names ---
# Here, the first two columns ("X" and "primary_gene") are not samples,
# so the sample names start from the 3rd column.
methylation_samples <- colnames(gene_level_methylation_csv)[-c(1,2)]
cat("\nGene-level methylation sample names (first few):\n")
print(head(methylation_samples))


# --- Find the Common Samples Among All Three Datasets ---
common_samples <- Reduce(intersect, list(short_transcriptomics, short_proteomics, methylation_samples))
cat("\nNumber of common samples among transcriptomics, proteomics, and methylation data:", 
    length(common_samples), "\n")
print(common_samples)
###############################################################
library(biomaRt)

# Connect to the Ensembl database for human genes.
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract peptide identifiers from the proteomics data.
# We assume the row names of proteomics_data are your peptide names.
peptides <- rownames(proteomics_data)
cat("Peptides in proteomics_data (first few):\n")
print(head(peptides))

# Query Ensembl using the "external_synonym" filter.
# This filter may retrieve gene mapping information if your peptide names appear
# as known synonyms/aliases. We request the official HGNC symbol.
mapping_results <- getBM(
  attributes = c("hgnc_symbol", "external_synonym"),
  filters = "external_synonym",
  values = peptides,
  mart = mart
)

# Check the mapping results.
cat("\nMapping results from biomaRt (first few rows):\n")
print(head(mapping_results))

# If the mapping worked, create a mapping table.
# Rename the column "external_synonym" to "peptide" so that we can merge by peptide.
if(nrow(mapping_results) > 0){
  mapping_table <- mapping_results
  colnames(mapping_table)[colnames(mapping_table) == "external_synonym"] <- "peptide"
} else {
  stop("No mapping results were found using biomaRt with the 'external_synonym' filter. 
       Consider checking your peptide names or using a custom mapping file.")
}

# Prepare the proteomics data for merging:
# Add the peptide names as a new column so they can be merged.
proteomics_data$peptide <- rownames(proteomics_data)

# Merge proteomics data with the mapping table.
proteomics_annotated <- merge(proteomics_data, mapping_table, by = "peptide", all.x = TRUE)
################

library(dplyr)

# Ensure the gene_name column is present. This will add it (or overwrite if needed)
proteomics_annotated <- proteomics_annotated %>%
  mutate(gene_name = hgnc_symbol)

# (Optional) Check that gene_name is now present:
print(colnames(proteomics_annotated))

# Aggregate duplicate entries by gene_name.
aggregated_proteomics <- proteomics_annotated %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

# Convert the result to a data frame.
aggregated_proteomics_df <- as.data.frame(aggregated_proteomics)

# Remove any rows with missing gene names
aggregated_proteomics_df <- aggregated_proteomics_df[!is.na(aggregated_proteomics_df$gene_name), ]

# Now set row names to gene_name.
rownames(aggregated_proteomics_df) <- aggregated_proteomics_df$gene_name

# (Optional) Remove the gene_name column if it's redundant.
 aggregated_proteomics_df$gene_name <- NULL

# Inspect the aggregated data.
cat("\nAggregated proteomics data (first few rows):\n")
head(aggregated_proteomics_df)
# Define the output file path. Note: Using forward slashes is preferred.
output_file <- "C:/Users/C11_LAB/Downloads/multiomics/aggregated_proteomics.csv"

# Save the data with row names (gene names)
write.csv(aggregated_proteomics_df, file = output_file, row.names = TRUE)

cat("Data successfully saved to", output_file, "\n")
