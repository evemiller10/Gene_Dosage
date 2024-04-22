convert_TCGAdata_to_Seg_CN <- function(project_name) {
  # Downloading necessary packages
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(dplyr)
  library(AnnotationDbi)
  library(doMC)
  
  # Source the functions for pulling in and converting TCGA data
  source ('/home/evem/CN_Seg/Scripts/Pull.convert.TCGA.functions.R')
  
  # Pulling in data
  TCGA.data <- pull_data_TCGA(project = project_name)
  
  #TCGA.data <- readRDS(file = "/home/evem/CN_Seg/Outputs/lgg.data.rds")
  
  # Save the data with the specific dataset name
  saveRDS(TCGA.data, file = paste0("/home/evem/CN_Seg/TCGA_Data/", project_name, ".data.rds"))
  
  TCGA.data$copynumber[!grepl("chr",TCGA.data$copynumber$Chromosome),] -> TCGA.data$copynumber
  table(factor(TCGA.data$copynumber$Chromosome))
  
  # Converting copy number TCGA.data from data frame to GRanges object
  cn.grange <- df_to_GRanges(TCGA.data$copynumber)
  
  # Manually adding CNcall metadata
  # >0.2= Amplification
  # <-0.2= Deletion
  # If other= Normal
  cn.grange$Call <- ifelse(mcols(cn.grange)$Segment_Mean > 0.2, 'Amplification',
                           ifelse(mcols(cn.grange)$Segment_Mean < -0.2, 'Deletion', 'Normal'))
  
  # Renaming Call metadata column as CNcall
  names(mcols(cn.grange))[names(mcols(cn.grange)) == 'Call'] <- 'CNcall'
  
  # Renaming Sample metadata column
  names(mcols(cn.grange))[names(mcols(cn.grange)) == 'sample.ids'] <- 'Sample'
  
  
  # Converts expression data from SummarizedExperiment object to GRanges object
  # Sample IDs not needed
  exp.grange <- rowRanges(TCGA.data$expression)
  
  
  # Source the function script for use below
  source ('/home/evem/CN_Seg/Scripts/exp2cn.source.R')            
  
  # creates a granges object with chromosomes binned at the sizes of 1Mb
  chromosome.map <- read.delim("./Archive/chromosome.map")
  
  #creates a granges object with certain custom regions of interest
  #keep.extra.columns = TRUE tells the function to keep all columns from the data frame that aren't used to construct the data frame
  #ignore.strand= TRUE tells the function to ignore the strand information when creating granges object
  select.features <- read.delim("./Archive/select.features.txt")
  makeGRangesFromDataFrame(chromosome.map, keep.extra.columns = TRUE,ignore.strand=TRUE) -> chr.grange
  makeGRangesFromDataFrame(select.features, keep.extra.columns = TRUE,ignore.strand=TRUE) -> select.grange
  
  # Renaming metadata columns within cn.grange so works within segMatrix function
  names(mcols(cn.grange)) <- sub("Segment_Mean", "Seg.CN", names(mcols(cn.grange)))
  names(mcols(cn.grange)) <- sub("Num_Probes", "Num.markers", names(mcols(cn.grange)))
  
  # install.packages("./diffloop_1.10.0.tar.gz", repos = NULL, type="source", dependencies = T)
  # Renaming seqlevels in exp.grange so they match cn.grange
  seqlevels(exp.grange)
  ## Rename 'seq2' to 'chr2' with a named vector.
  exp.grange.2 <- renameSeqlevels(exp.grange, sub("chr", "", seqlevels(exp.grange)))
  # Rename seqlevels in chr.grange and select.grange
  chr.grange.2 <- renameSeqlevels(chr.grange, sub("chr", "", seqlevels(chr.grange)))
  select.grange.2 <- renameSeqlevels(select.grange, sub("chr", "", seqlevels(select.grange)))
  
  #generates a list of matrices where individuals represented by columns and genes as rows
  #contains 4 slots:
  #seg.means= average CN ration within each segment
  #CN
  #CNcall= categorical calls representing aplifications/deletions etc
  #num.mark= number of markers within each segment
  segMatrix(exp.grange.2, cn.grange) -> seg.list
  segMatrix(chr.grange.2, cn.grange) -> seg.chr.list
  segMatrix(select.grange.2, cn.grange) -> seg.select.list
  
  # Remove NA values from each element of the list
  seg.select.list$seg.means <- na.omit(seg.select.list$seg.means)
  seg.select.list$CNcall <- na.omit(seg.select.list$CNcall)
  seg.select.list$num.mark <- na.omit(seg.select.list$num.mark)
  
  # Extract assay data from TCGA.data$expression 
  exp.assay.data <- assay(TCGA.data$expression, "tpm_unstrand")
  
  
  #create an index select and match genes
  #matches the row names of seg.list$CNcall and TCGA expression data
  #subsets data to only keep rows where the corresponding row in index is not missing
  ### Think this is not correctly matching Ensembl gene IDs ###
  #match (rownames(seg.list$CNcall), rownames(TCGA.data$expression)) -> index
  ### Trying an alternative method to match gene IDs ###
  index <- match(rownames(seg.list$CNcall), rownames(exp.assay.data))
  exp.grange.2[which(!is.na(index))] -> exp.grange.2
  exp.assay.data[index[!is.na(index)],] -> exp.assay.data
  seg.list$CNcall[which(!is.na(index)),] -> seg.list$CNcall
  seg.list$seg.means[which(!is.na(index)),] -> seg.list$seg.means
  
  # Correctly matches column names
  # TCGA data had longer sample IDs within CN data
  # Removes additional information from sample IDs to match expression data
  make.unique(sub("^([^-]*-[^-]*-[^-]*).*", "\\1", colnames(seg.list$CNcall))) -> colnames(seg.list$CNcall)
  make.unique(sub("^([^-]*-[^-]*-[^-]*).*", "\\1", colnames(seg.list$seg.means))) -> colnames(seg.list$seg.means)
  make.unique(sub("^([^-]*-[^-]*-[^-]*).*", "\\1", colnames(exp.assay.data))) -> colnames(exp.assay.data)
  match(colnames(seg.list$CNcall), colnames(exp.assay.data)) -> index
  
  #match the copy number sample IDs to the available expression data column IDs
  seg.list$CNcall[,!is.na(index)] -> matched.cn
  # seg.list$seg.means[,which(!is.na(index))] -> matched.seg.means
  seg.list$seg.means[,!is.na(index)] -> matched.seg.means
  exp.assay.data[,index[!is.na(index)]] -> matched.exp.assay.data
  
  # Checks column names are identical
  identical(colnames(matched.exp.assay.data), colnames(matched.seg.means))
  identical(colnames(matched.exp.assay.data), colnames(matched.cn))
  
  # Checks row names are identical
  identical(rownames(matched.exp.assay.data), rownames(matched.seg.means))
  identical(rownames(matched.exp.assay.data), rownames(matched.cn))
  
  # Filtering expression genes/copynumber
  dim(matched.exp.assay.data)
  apply(matched.exp.assay.data, 1, gp.style.filter, fold.change=3, delta=10, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE) -> index.genes
  filt.matched.exp.assay.data <- matched.exp.assay.data[index.genes,]
  filt.matched.seg.means <- matched.seg.means[index.genes,]
  filt.matched.cn <- matched.cn[index.genes,]
  
  dim(filt.matched.exp.assay.data)
  dim(filt.matched.seg.means)
  dim(filt.matched.cn)
  
  # set the minumum number of samples possessing copy number changes 
  min.number.samples <- 5
  
  saveRDS(matched.cn, file = "/home/evem/CN_Seg/TCGA_Outputs/matched.cn.rds")
  saveRDS(matched.seg.means, file = "/home/evem/CN_Seg/TCGA_Outputs/matched.seg.mean.rds")
  saveRDS(matched.exp.assay.data, file = "/home/evem/CN_Seg/TCGA_Outputs/matched.exp.assay.data.rds")
  
  saveRDS(filt.matched.cn, file = "/home/evem/CN_Seg/TCGA_Outputs/filt.matched.cn.rds")
  saveRDS(filt.matched.seg.means, file = "/home/evem/CN_Seg/TCGA_Outputs/filt.matched.seg.mean.rds")
  saveRDS(filt.matched.exp.assay.data, file = "/home/evem/CN_Seg/TCGA_Outputs/filt.matched.exp.assay.data.rds")
  
  # set up normal margins
  # calculate normal values
  ### This is where the error occurs ###
  matched.cn <- readRDS(file = "/home/evem/CN_Seg/TCGA_Outputs/matched.cn.rds")
  matched.exp.assay.data <- readRDS(file = "/home/evem/CN_Seg/TCGA_Outputs/matched.exp.assay.data.rds")
  normal.values <- calculate.normal.vals(matched.cn, matched.exp.assay.data, cores = 2)
  
  # Create variables with project abbreviation in their names
  # Adds each variable to the current environment
  assign(paste0("cn.", project_name, ".grange"), cn.grange, envir = .GlobalEnv)
  assign(paste0("exp.", project_name, ".assay.data"), exp.assay.data, envir = .GlobalEnv)
  assign(paste0("exp.", project_name, ".grange"), exp.grange, envir = .GlobalEnv)
  assign(paste0("exp.", project_name, ".grange.2"), exp.grange.2, envir = .GlobalEnv)
  assign(paste0("filt.matched.", project_name, ".cn"), filt.matched.cn, envir = .GlobalEnv)
  assign(paste0("filt.matched.", project_name, ".exp.assay.data"), filt.matched.exp.assay.data, envir = .GlobalEnv)
  assign(paste0("filt.matched.", project_name, ".seg.means"), filt.matched.seg.means, envir = .GlobalEnv)
  assign(paste0("matched.", project_name, ".exp.assay.data"), matched.exp.assay.data, envir = .GlobalEnv)
  assign(paste0("matched.", project_name, ".seg.means"), matched.seg.means, envir = .GlobalEnv)
  assign(paste0(project_name, ".normal.values"), normal.values, envir = .GlobalEnv)
  assign(paste0("seg.list.", project_name), seg.list, envir = .GlobalEnv)
  assign(paste0("select.", project_name, ".features"), select.features, envir = .GlobalEnv)
  assign(paste0("select.", project_name, ".grange"), select.grange, envir = .GlobalEnv)
  assign(paste0("select.", project_name, ".grange.2"), select.grange.2, envir = .GlobalEnv)
}

convert_TCGAdata_to_Seg_CN("TCGA-LGG")
convert_TCGAdata_to_Seg_CN("TCGA-ACC")
convert_TCGAdata_to_Seg_CN("TCGA-SARC")
convert_TCGAdata_to_Seg_CN("TCGA-KIRP")