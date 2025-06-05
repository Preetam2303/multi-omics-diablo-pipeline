install.packages("BiocManager")
BiocManager::install(c("minfi", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "IlluminaHumanMethylation450kmanifest"))
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
BiocManager::install("sesame")
library(sesame)

data_dir <- "C:/Users/C11_LAB/Downloads/gdc_downloads"
beta_values <- openSesame(data_dir)

# Replace this with your actual root directory
root_dir <- "C:/Users/C11_LAB/Downloads/gdc_downloads"

# List all relevant methylation beta files recursively
files <- list.files(
  path = root_dir, 
  pattern = "\\.methylation_array\\.sesame\\.level3betas$", 
  recursive = TRUE, 
  full.names = TRUE
)

# Check if the list is correct:
print(files)

# Check if the directory exists
dir.exists(root_dir)

# List all files in the directory (including subdirectories)
all_files <- list.files(path = root_dir, recursive = TRUE, full.names = TRUE)
print(all_files)
#################
# Load necessary libraries
library(data.table)
library(dplyr)

# 1. Define the root directory (adjust to your system path)
root_dir <- "C:/Users/C11_LAB/Downloads/gdc_downloads"

# 2. List all files recursively from the root directory
all_files <- list.files(
  path = root_dir, 
  recursive = TRUE, 
  full.names = TRUE
)

# 3. Filter to include only the methylation beta files while excluding files in "logs" directories.
#    Adjust the regular expression if your file naming is slightly different.
files_filtered <- all_files[
  grepl("methylation_array\\.sesame\\.level3betas\\.txt$", all_files, ignore.case = TRUE) & 
    !grepl("/logs/", all_files, ignore.case = TRUE)
]

# Check the filtered file list
print(files_filtered)

# 4. Define a function to read each methylation file and assign its sample name.
#    In this case, we assume the sample ID is the name of the folder containing the file.
read_methyl_file <- function(file) {
  # Extract the sample name from the parent folder name
  sample_id <- basename(dirname(file))
  
  # Read the file: assumes two columns (CpG_ID and beta value)
  dt <- fread(file, header = FALSE, na.strings = "N/A")
  
  # Rename the columns: first column as "CpG_ID" and second column as the sample_id
  setnames(dt, c("CpG_ID", sample_id))
  
  return(dt)
}

# 5. Read each filtered file using lapply to get a list of data tables
methyl_list <- lapply(files_filtered, read_methyl_file)

# 6. Merge all data tables on the common "CpG_ID" column.
#    full_join ensures every CpG site is kept even if data is missing in some samples.
combined_df <- Reduce(function(x, y) full_join(x, y, by = "CpG_ID"), methyl_list)

# Optionally check the combined data frame
head(combined_df)

# 7. Convert the merged data frame to a matrix for downstream analysis.
#    Rows will be the CpG probe IDs and columns correspond to samples.
beta_matrix <- as.matrix(combined_df[,-1])
rownames(beta_matrix) <- combined_df$CpG_ID

# Check the resulting matrix dimensions
dim(beta_matrix)
# Check the structure of the combined data frame
str(combined_df)

####################
# Convert all columns except the first to numeric
combined_df_numeric <- combined_df
combined_df_numeric[,-1] <- lapply(combined_df_numeric[,-1], as.numeric)

# Verify the conversion by checking the structure again
str(combined_df_numeric)
# Convert columns 2 through ncol(combined_df_numeric) to numeric
combined_df_numeric[, 2:ncol(combined_df_numeric)] <- 
  lapply(combined_df_numeric[, 2:ncol(combined_df_numeric), with = FALSE], as.numeric)
##
# Define a helper function to list non-numeric entries in a vector
find_non_numeric <- function(x) {
  # Convert to numeric but capture the original values that yield NA (ignoring those that are "NA" strings if you expect them)
  non_numeric <- x[is.na(as.numeric(x)) & !is.na(x) & x != "NA"]
  unique(non_numeric)
}

# Check the first column (excluding the CpG_ID) for non-numeric entries
problematic_values <- find_non_numeric(combined_df_numeric[[2]])
print(problematic_values)
##
# Check columns 3 to 5 for non-numeric entries
problematic_values_list <- lapply(3:5, function(i) {
  list(
    column = colnames(combined_df_numeric)[i],
    non_numeric = find_non_numeric(combined_df_numeric[[i]])
  )
})
print(problematic_values_list)
na_counts <- combined_df_numeric[, lapply(.SD, function(col) sum(is.na(col))), .SDcols = 2:ncol(combined_df_numeric)]
print(na_counts)
library(ggplot2)
str(beta_df)
beta_df <- data.frame(lapply(beta_df, function(x) as.numeric(as.character(x))))
str(beta_df)  # Verify that all columns are now numeric
boxplot(beta_df, outline = FALSE, main = "Beta Value Distribution (First 10 Samples)")
dim(beta_matrix)
head(rownames(beta_matrix))
head(colnames(beta_matrix))
# Save just the beta_matrix
saveRDS(beta_matrix, file = "beta_matrix.rds")
beta_matrix <- readRDS("beta_matrix.rds")
#########################################
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann)
# Convert beta_matrix to a data frame and create a CpG_ID column
beta_df <- as.data.frame(beta_matrix)
beta_df$CpG_ID <- rownames(beta_matrix)
# Ensure the annotation data frame has a CpG_ID column based on its row names
ann_df <- as.data.frame(ann)
ann_df$CpG_ID <- rownames(ann_df)
##
# Merge beta_df with the annotation information using the CpG_ID as key
merged_data <- merge(beta_df, ann_df[, c("CpG_ID", "UCSC_RefGene_Name")], 
                     by = "CpG_ID", all.x = TRUE)

#####################################
# Extract the primary gene for each CpG
merged_data$primary_gene <- sapply(merged_data$UCSC_RefGene_Name, function(x) {
  if (!is.na(x) && nchar(x) > 0) {
    strsplit(x, ";")[[1]][1]  # Taking the first gene listed
  } else {
    NA
  }
})
#################
library(dplyr)
str(merged_data)
# Convert columns 2 to 512 to numeric (adjust indices if needed)
merged_data[, 2:512] <- lapply(merged_data[, 2:512], function(x) as.numeric(x))

# Check the structure again to ensure they are numeric
str(merged_data)
library(dplyr)

gene_level_methylation <- merged_data %>%
  group_by(primary_gene) %>%
  summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

# Display the first few rows of the gene-level summary
head(gene_level_methylation)
# See the names of the aggregated numeric columns (i.e., the sample IDs)
colnames(gene_level_methylation)

# Get a summary of the first few numeric columns
summary(gene_level_methylation[, 2:6])

saveRDS(gene_level_methylation, file = "gene_level_methylation.rds")

# Later, load it with:
gene_level_methylation <- readRDS("gene_level_methylation.rds")
#############################
# Load the sample sheet from the specified path
sample_sheet <- read.delim("C:/Users/C11_LAB/Downloads/gdc_sample_sheet.2025-04-14 (1).tsv",
                           header = TRUE, stringsAsFactors = FALSE)

# Inspect the sample sheet structure and first few rows
str(sample_sheet)
head(sample_sheet)

# Assuming the columns in the sample sheet are named "File.ID" and "Sample.ID"
mapping <- sample_sheet$Sample.ID
names(mapping) <- sample_sheet$File.ID

# Check the mapping
head(mapping)

# Suppose beta_matrix currently has column names corresponding to the File.IDs
old_names <- colnames(beta_matrix)

# Look up the new names from your mapping vector:
new_names <- mapping[old_names]

# Replace the column names: if a file ID isn't found in the mapping (NA), keep the original name.
colnames(beta_matrix) <- ifelse(is.na(new_names), old_names, new_names)

# Verify the updated column names
head(colnames(beta_matrix))
missing_names <- old_names[is.na(new_names)]
print(missing_names)
###################
# Extract the current column names, excluding the first "primary_gene" column
old_gene_names <- colnames(gene_level_methylation)[-1]

# Look up the new names using the mapping vector
new_gene_names <- mapping[old_gene_names]

# Replace file-based column names with TCGA sample IDs (keeping original names if mapping is not found)
colnames(gene_level_methylation)[-1] <- ifelse(is.na(new_gene_names), old_gene_names, new_gene_names)

# Check the updated column names
head(colnames(gene_level_methylation))

################
# Check the structure of the data
str(gene_level_methylation)

# View the first few rows
head(gene_level_methylation)

# Optionally, get a summary of the numeric columns (excluding the 'primary_gene' column)
summary(gene_level_methylation[, -1])
#############
# Save the data as an RDS file
saveRDS(gene_level_methylation, file = "gene_level_methylation.rds")

# Later, load it with:
# gene_level_methylation <- readRDS("gene_level_methylation.rds")

###################
# Load the RPPA data from the CSV file
rppa <- read.csv("C:/Users/C11_LAB/Downloads/TCGA_HNSC_RPPA.csv",
                 header = TRUE, stringsAsFactors = FALSE)

# Inspect the first few rows and structure if needed
str(rppa)
head(rppa)

# In the RPPA data, the first six columns are metadata (e.g., X, AGID, lab_id, etc.)
# and the sample data start from the 7th column onward.
# Extract the sample IDs from these columns:
rppa_sample_ids_raw <- colnames(rppa)[7:ncol(rppa)]

# The RPPA sample IDs are in a period-separated format (e.g., "TCGA.CV.7423.01A").
# Convert them into hyphen-separated IDs to match your gene_level_methylation sample IDs
# (e.g., "TCGA-CV-7423-01A").
rppa_sample_ids <- gsub("\\.", "-", rppa_sample_ids_raw)

# For your methylation data (gene_level_methylation), the first column is 'primary_gene'
# and the rest are sample IDs that are already hyphenated.
methyl_sample_ids <- colnames(gene_level_methylation)[-1]

# Find the common sample IDs between the two datasets.
common_samples <- intersect(rppa_sample_ids, methyl_sample_ids)

# Count the number of common samples.
common_count <- length(common_samples)
cat("Number of common sample IDs:", common_count, "\n")

# Optionally, print the list of common sample IDs.
print(common_samples)
