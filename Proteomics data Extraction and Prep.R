# Define the path to your tar.gz file. Notice the forward slashes.
tar_file <- "C:/Users/C11_LAB/Downloads/gdc_download_20250421_055731.977352.tar.gz"

# Define an extraction directory. You can create one, or use an existing folder.
extraction_dir <- "C:/Users/C11_LAB/Downloads/extracted"

# Create the extraction directory if it doesn't exist
if (!dir.exists(extraction_dir)) {
  dir.create(extraction_dir)
}

# Extract the archive
untar(tar_file, exdir = extraction_dir)
# List all files in the extraction directory to see what's inside
extracted_files <- list.files(extraction_dir, recursive = TRUE, full.names = TRUE)
print(extracted_files)
###############
# Define the base extraction directory
extraction_dir <- "C:/Users/C11_LAB/Downloads/extracted"

# List only TSV files matching the RPPA_data naming convention; this excludes MANIFEST.txt automatically.
tsv_files <- list.files(
  path = extraction_dir, 
  pattern = "_RPPA_data\\.tsv$", 
  recursive = TRUE, 
  full.names = TRUE
)

# Check the list to confirm that MANIFEST.txt is not included
print(tsv_files)
##################
# Use the first file in the list as an example
example_file <- tsv_files[1]
# Read the first TSV file; adjust header and separator if necessary
example_data <- read.delim(example_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Check the structure of the data frame
str(example_data)

# See the first few rows of the data frame
head(example_data)

# Optionally, check the dimensions to see how many genes/proteins are reported
dim(example_data)
#################################################
# Use the previously defined vector `tsv_files` of your file paths.

data_list <- lapply(tsv_files, function(file) {
  # Read the TSV file
  df <- read.delim(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract the gene/protein identifiers and expression values.
  # Here we use 'peptide_target' as the gene identifier,
  # and 'protein_expression' as the expression measurement.
  targets <- df$peptide_target
  expr_values <- df$protein_expression
  
  # Create a new data frame with row names set to the peptide targets;
  # the column is the expression values.
  sample_df <- data.frame(Expression = expr_values, row.names = targets, stringsAsFactors = FALSE)
  
  # Derive a sample ID from the file name by removing the suffix.
  sample_id <- gsub("_RPPA_data\\.tsv$", "", basename(file))
  colnames(sample_df) <- sample_id
  
  return(sample_df)
})

#######################
# Convert row names (gene identifiers) into a column called 'peptide_target'
data_list_mod <- lapply(data_list, function(df) {
  df_mod <- cbind(peptide_target = rownames(df), df)
  return(df_mod)
})

# Merge all data frames by the 'peptide_target' column.
# The 'all = TRUE' argument makes sure that even if a gene is missing in one sample,
# it will appear in the final merged data with an NA.
combined_data <- Reduce(function(x, y) merge(x, y, by = "peptide_target", all = TRUE), data_list_mod)

# Reassign row names using the 'peptide_target' and remove the now redundant column.
rownames(combined_data) <- combined_data$peptide_target
combined_data$peptide_target <- NULL

# Inspect the combined data frame to ensure genes are rows and sample identifiers are columns.
head(combined_data)
##############
output_file <- "C:/Users/C11_LAB/Downloads/proteomics_combined.csv"
write.csv(combined_data, file = output_file, row.names = TRUE)
cat("Data successfully written to", output_file, "\n")
