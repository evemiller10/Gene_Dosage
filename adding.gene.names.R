add_gene_names <- function(df) {
  # Downloading package for gene name conversion
  library(biomaRt)
  
  # Getting ensembl IDs from package biomaRt
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
  
  # Removing version number from ensembl IDs in exp assay data df
  rownames(df) <- sub("\\.[0-9]+$", "", rownames(df))
  
  # Getting gene names from biomaRt
  gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                      filters = "ensembl_gene_id", 
                      values = rownames(df), 
                      mart = ensembl)
  
  # Convert row names to a column
  df$ensembl_gene_id <- rownames(df)
  
  # Add original rownames as a new column
  df$original_row_names <- rownames(df)
  df$row_number <- 1:nrow(df)
  
  # Merge df with gene_names based on ensembl_gene_id
  merged_data <- merge(df, gene_names, by = "ensembl_gene_id", all.x = TRUE)
  
  # Remove the temporary ensembl_gene_id column
  merged_data$ensembl_gene_id <- NULL
  
  # Rename the external_gene_name column to GeneName
  colnames(merged_data)[which(names(merged_data) == "external_gene_name")] <- "GeneName"
  
  # Order by row number to restore the original order of the rows
  merged_data <- merged_data[order(merged_data$row_number), ]
  
  # Set the original row names as row names of merged data
  rownames(merged_data) <- merged_data$original_row_names
  
  # Remove the temporary original_row_names column
  merged_data$original_row_names <- NULL
  merged_data$row_number <- NULL
  
  # Return the updated data frame
  return(merged_data)
}
