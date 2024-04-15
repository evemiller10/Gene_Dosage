# Downloading necessary packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(GenomicRanges)
library(GenomicFeatures)
library(dplyr)
library(AnnotationDbi)
library(doMC)

# Singular function to pull in all necessary TCGA data
pull_data_TCGA <- function(
    project = "",
    workflow.type = "STAR - Counts",
    save.filename = ""
){
  
  # Downloading Expression Data
  query.exp <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    #below line was producing error as unused argument> don't think it is needed here
    #legacy = FALSE,
    data.type = "Gene Expression Quantification",
    workflow.type = workflow.type)
  GDCdownload(query.exp)
  expression <- GDCprepare(query = query.exp)
  
  
  # Downloading Copy Number Data
  query.cn <- GDCquery(
    project = project,
    data.category = "Copy Number Variation",
    data.type = "Copy Number Segment")
  GDCdownload(query.cn)
  copynumber <- GDCprepare(query.cn)
  
  # Downloading Mutation Data
  query.mutation <- GDCquery(
    project = project,
    data.category = "Simple Nucleotide Variation",
    access = "open",
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  GDCdownload(query.mutation)
  mutational <- GDCprepare(query.mutation)
  
  # Downloading Clinical Data
  query.clinical <- GDCquery(
    project = project,
    #file.type = "xml",
    data.category = "Clinical",
    #legacy = FALSE,
    data.type = "Clinical Supplement",
    data.format = "bcr xml")
  GDCdownload(query.clinical)
  clinical <- GDCprepare_clinic(query = query.clinical, clinical.info = "patient")
  
  # Returns the elements as a list
  return(list(expression = expression,
              clinical = clinical,
              copynumber = copynumber,
              mutational = mutational))
}

# Pulling in GBM data
TCGA.GBMdata <- pull_data_TCGA(project = "TCGA-GBM", save.filename = "/home/evem/CN_Seg/Outputs")

# Function to convert dataframe to Granges 
# Used for CN and mutational data
df_to_GRanges <- function(df = "", seqnames.field = "Chromosome", start.field = c("Start", "Start_Position"), end.field = c("End", "End_Position"), strand.field = "Strand") {
  granges <- makeGRangesFromDataFrame(df,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqinfo = NULL,
                                      seqnames.field = seqnames.field,
                                      start.field = start.field,
                                      end.field = end.field,
                                      strand.field = strand.field)
  
}
cn.gbm.grange <- df_to_GRanges(df = TCGA.GBMdata$copynumber)

# Manually adding CNcall metadata
# >0.2= Amplification
# <-0.2= Deletion
# If other= Normal
cn.gbm.grange$Call <- ifelse(mcols(cn.gbm.grange)$Segment_Mean > 0.2, 'Amplification',
                             ifelse(mcols(cn.gbm.grange)$Segment_Mean < -0.2, 'Deletion', 'Normal'))

# Converts expression data from SummarizedExperiment object to GRanges object
exp.gbm.grange <- rowRanges(TCGA.GBMdata$expression)

# Source the function script for use below
source ('./Scripts/exp2cn.source.R')

# creates a granges object with chromosomes binned at the sizes of 1Mb
chromosome.map <- read.delim("./Archive/chromosome.map")

#creates a granges object with certain custom regions of interest
#keep.extra.columns = TRUE tells the function to keep all columns from the data frame that aren't used to construct the data frame
#ignore.strand= TRUE tells the function to ignore the strand information when creating granges object
select.features <- read.delim("./Archive/select.features.txt")
makeGRangesFromDataFrame(chromosome.map, keep.extra.columns = TRUE,ignore.strand=TRUE) -> chr.grange
makeGRangesFromDataFrame(select.features, keep.extra.columns = TRUE,ignore.strand=TRUE) -> select.grange

# Renaming metadata columns within cn.gbm.grange so works within segMatrix function
names(mcols(cn.gbm.grange)) <- sub("Segment_Mean", "Seg.CN", names(mcols(cn.gbm.grange)))
names(mcols(cn.gbm.grange)) <- sub("Num_Probes", "Num.markers", names(mcols(cn.gbm.grange)))

#generates a list of matrices where individuals represented by columns and genes as rows
#contains 4 slots:
#seg.means= average CN ration within each segment
#CN
#CNcall= categorical calls representing aplifications/deletions etc
#num.mark= number of markers within each segment
### No error but not sure it is creating seg.list.gbm correctly ###
segMatrix(exp.gbm.grange, cn.gbm.grange) -> seg.list.gbm
segMatrix(chr.grange,cn.gbm.grange) -> seg.chr.list
segMatrix(select.grange,cn.gbm.grange) -> seg.select.list

# Remove NA values from each element of the list
seg.select.list$seg.means <- na.omit(seg.select.list$seg.means)
seg.select.list$CNcall <- na.omit(seg.select.list$CNcall)
seg.select.list$num.mark <- na.omit(seg.select.list$num.mark)

#create an index select and match genes
#matches the row names of seg.list.gbm$CNcall and TCGA expression data
#subsets data to only keep rows where the corresponding row in index is not missing
### Think this is not correctly matching Ensembl gene IDs ###
#match (rownames(seg.list.gbm$CNcall), rownames(TCGA.GBMdata$expression)) -> index
### Trying an alternative method to match gene IDs ###
index <- match(rownames(seg.list.gbm$CNcall), rownames(TCGA.GBMdata$expression))
exp.gbm.grange[which(!is.na(index))] -> exp.gbm.grange
TCGA.GBMdata$expression[index[!is.na(index)],] -> TCGA.GBMdata$expression
seg.list.gbm$CNcall[which(!is.na(index)),] -> seg.list.gbm$CNcall
seg.list.gbm$seg.means[which(!is.na(index)),] -> seg.list.gbm$seg.means

#match the expression data column IDs to the available copy number sample IDs
### match(colnames(seg.list.gbm$CNcall), colnames(TCGA.GBMdata$expression)) -> index ###
index <- match(colnames(seg.list.gbm$CNcall), colnames(TCGA.GBMdata$expression))
TCGA.GBMdata$expression[,index[!is.na(index)]] -> select.exp.data

#match the copy number sample IDs to the available expression data column IDs
seg.list.gbm$CNcall -> matched.cn
seg.list.gbm$seg.means -> matched.seg.means


# set the minumum number of samples possessing copy number changes 
min.number.samples <- 5

# set up normal margins
# calculate normal values
### This is where the error occurs ###
normal.values <- calculate.normal.vals(matched.cn, select.exp.data, cores =20)
sample.number <- 10