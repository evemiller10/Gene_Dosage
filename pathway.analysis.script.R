#Loading in necessary packages 
library(fgsea)
library(data.table)
library(biomaRt)
library(ggplot2)

#Function for creating coefficient dataframe for rank conversion
rank_LM_results <- function(tumour_type_LM_results_scaled, tumour_type) {
  # Convert to data frame
  tumour_type_LM_results_scaled.df <- as.data.frame(tumour_type_LM_results_scaled)
  
  # Extract lm.seg.means.coef column
  coefs <- tumour_type_LM_results_scaled.df$lm.seg.means.coef
  
  # Order coefficients from highest to lowest, preserving NA values
  ordered_coefs <- sort(coefs, decreasing = TRUE, na.last = TRUE)
  
  # Create a matrix with ordered coefficients
  expression_matrix <- matrix(ordered_coefs, ncol = 1)
  
  # Set row names to Ensembl IDs from the original data frame
  rownames(expression_matrix) <- rownames(tumour_type_LM_results_scaled.df)
  
  # Optionally, assign column names
  colnames(expression_matrix) <- "lm.seg.means.coef"
  
  # Convert matrix to data frame
  expression_df <- as.data.frame(expression_matrix)
  
  # Add gene names (assuming add_gene_names function is defined elsewhere)
  expression_df <- add_gene_names(expression_df)
  
  # Assign the data frame to a variable name that includes the tumour type
  assign(paste0("expression_df_", tumour_type), expression_df, envir = .GlobalEnv)
  
  return(expression_df)
}

expression_df_ACC <- rank_LM_results(ACC_LM_results_scaled.df, "ACC")
expression_df_ESCA <- rank_LM_results(ESCA_LM_results_scaled.df, "ESCA")
expression_df_LGG <- rank_LM_results(LGG_LM_results_scaled.df, "LGG")
expression_df_GBM <- rank_LM_results(GBM_LM_results_scaled.df, "GBM")
expression_df_BLCA <- rank_LM_results(BLCA_LM_results_scaled.df, "BLCA")
expression_df_KIRP <- rank_LM_results(KIRP_LM_results_scaled.df, "KIRP")
expression_df_LIHC <- rank_LM_results(LIHC_LM_results_scaled.df, "LIHC")
expression_df_SARC <- rank_LM_results(SARC_LM_results_scaled.df, "SARC")
expression_df_PAAD <- rank_LM_results(PAAD_LM_results_scaled.df, "PAAD")
expression_df_TGCT <- rank_LM_results(TGCT_LM_results_scaled.df, "TGCT")


################################################################################


#Downloading pathway data sets from MSigDB
#c5 data
gmt_url <-"https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.all.v2023.2.Hs.symbols.gmt"

#Hallmark genes data
gmt_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt"

#c6 Oncogenic Gene Set
gtm_url <-  "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c6.all.v2023.2.Hs.symbols.gmt"

#Setting directory to save the file in
destination_dir <- "/home/evem/CN_Seg/H.all.v2023.2.Hs.symbols.gmt"
destination_dir <- "home/evem/CN_Seg"
destination_file <- "c6.all.v2023.2.Hs.symbols.gmt"
destfile <- file.path(destination_dir, destination_file)
download.file("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c6.all.v2023.2.Hs.symbols.gmt", destfile)

# Download the file
download.file(gmt_url, destfile = destination_dir, method = "auto")

#Converting gmt file to format needed for pathway analysis
library(qusage)

# Define the path to the downloaded GMT file
gmt_file <- "/home/evem/CN_Seg/home/evem/CN_Seg/c6.all.v2023.2.Hs.symbols.gmt"
gmt_file <- "/home/evem/CN_Seg/H.all.v2023.2.Hs.symbols.gmt"

# Read the contents of the GMT file
pathways <- read.gmt(gmt_file)

#Ensure pathways file is in correct format after being read
head(pathways)


################################################################################

# rankings is called temp.df

#Function for 
rank_expression_df <- function(expression_df, tumor_type) {
  # Step 1: Extract numeric values and assign names
  ranks <- as.numeric(expression_df[, -2])
  names(ranks) <- expression_df[, 2]
  
  # Step 2: Filter out non-finite values (NA, Inf, -Inf)
  ranks.filt <- ranks[is.finite(ranks) & !is.infinite(ranks)]
  
  # Step 3: Ensure unique names
  ranks.filt <- ranks.filt[!duplicated(names(ranks.filt))]
  
  # Dynamically create object name with tumor_type
  rank_filt_name <- paste("rank.filt", tumor_type, sep = "_")
  
  # Assign the filtered ranks.filt object to the dynamically named variable
  assign(rank_filt_name, ranks.filt, envir = .GlobalEnv)
  
  # Optionally, return the object name (useful for checking or further manipulation)
  return(rank_filt_name)
}

rank_expression_df(expression_df_ACC, "ACC")
rank_expression_df(expression_df_ESCA, "ESCA")
rank_expression_df(expression_df_LGG, "LGG")
rank_expression_df(expression_df_GBM, "GBM")
rank_expression_df(expression_df_BLCA, "BLCA")
rank_expression_df(expression_df_KIRP, "KIRP")
rank_expression_df(expression_df_LIHC, "LIHC")
rank_expression_df(expression_df_SARC, "SARC")
rank_expression_df(expression_df_PAAD, "PAAD")
rank_expression_df(expression_df_TGCT, "TGCT")

################################################################################

fgseaRes_ACC <- fgsea(pathways = pathways, 
                  stats    = rank.filt_ACC,
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes_ACC[order(padj), ],30)
fgseaRes_GBM_df <- as.data.frame(fgseaRes_GBM)
# Filter fgseaRes_ACC_df for the HALLMARK_E2F_TARGETS pathway
kras_targets_df <- fgseaRes_SARC_df[fgseaRes_SARC_df$pathway == "HALLMARK_KRAS_SIGNALING_UP", ]

# Extract leading edge genes for the HALLMARK_E2F_TARGETS pathway
leading_edge_genes_kras <- lapply(kras_targets_df$leadingEdge, function(x) unlist(strsplit(x, "/")))

# Flatten leading_edge_genes_e2f into a single vector
leading_edge_genes_kras_flat <- unlist(leading_edge_genes_kras)

# Convert gene names to a consistent case (e.g., uppercase) for both lists
leading_edge_genes_kras_flat <- toupper(leading_edge_genes_kras_flat)
# Assuming chrom1p_gene_names_upper is already defined from previous steps

# Check for any leading edge genes in chrom1p_gene_names for HALLMARK_E2F_TARGETS
matching_genes_kras_1p_SARC <- leading_edge_genes_kras_flat[leading_edge_genes_kras_flat %in% chrom1p_gene_names_upper]

fgseaRes_BLCA <- fgsea(pathways = pathways, 
                      stats    = rank.filt_BLCA,
                      minSize  = 15,
                      maxSize  = 500)

head(fgseaRes_BLCA[order(padj), ],30)

fgseaRes_ESCA <- fgsea(pathways = pathways, 
                      stats    = rank.filt_ESCA,
                      minSize  = 15,
                      maxSize  = 500)

head(fgseaRes_ESCA[order(padj), ],30)

fgseaRes_GBM <- fgsea(pathways = pathways, 
                      stats    = rank.filt_GBM,
                      minSize  = 15,
                      maxSize  = 500)

head(fgseaRes_GBM[order(padj), ],30)

fgseaRes_KIRP <- fgsea(pathways = pathways, 
                      stats    = rank.filt_KIRP,
                      minSize  = 15,
                      maxSize  = 500)

head(fgseaRes_KIRP[order(padj), ],30)

fgseaRes_LGG <- fgsea(pathways = pathways, 
                      stats    = rank.filt_LGG,
                      minSize  = 15,
                      maxSize  = 500)

head(fgseaRes_LGG[order(padj), ],30)

fgseaRes_LIHC <- fgsea(pathways = pathways, 
                      stats    = rank.filt_LIHC,
                      minSize  = 15,
                      maxSize  = 500)

head(fgseaRes_LIHC[order(padj), ],30)

fgseaRes_PAAD <- fgsea(pathways = pathways, 
                      stats    = rank.filt_PAAD,
                      minSize  = 15,
                      maxSize  = 500)

head(fgseaRes_PAAD[order(padj), ],30)

fgseaRes_SARC <- fgsea(pathways = pathways, 
                      stats    = rank.filt_SARC,
                      minSize  = 15,
                      maxSize  = 500)

head(fgseaRes_SARC[order(padj), ],30)

fgseaRes_TGCT <- fgsea(pathways = pathways, 
                      stats    = rank.filt_TGCT,
                      minSize  = 15,
                      maxSize  = 500)

head(fgseaRes_TGCT[order(padj), ],30)
