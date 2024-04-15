### Coverts TCGA expression, methylation, copy number, and mutation data into GRanges objects ###

# Required Packages
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(sesame)
library(sesameData)
library(SummarizedExperiment)
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationDbi)

# Pulling TCGA CN data and converting to a GRanges object
pull.cn.data.TCGA <- function(project = "", workflow.type = "STAR - Counts", save.filename = ""){
  query.cn <- GDCquery(
    project = project,
    data.category = "Copy Number Variation",
    data.type = "Copy Number Segment")
  GDCdownload(query.cn)
  copynumber.data <- GDCprepare(query.cn)
  return(copynumber.data)
}

# Manually adding CNcall metadata
cn_grange.gbm$CNcall <- ifelse(mcols(cn_grange.gbm)$Segment_Mean > 0.2, 'Amplification',
                               ifelse(mcols(cn_grange.gbm)$Segment_Mean < -0.2, 'Deletion', 'Normal'))

# Converting to GRanges object from DataFrame  
makeGRangesFromDataFrame(cn.gbm,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqinfo = NULL,
      seqnames.field = "Chromosome",
      start.field = c("Start", "Start_Position"),
      end.field = c("End", "End_Position"),
      strand.field = "Strand") -> cn_grange.gbm

# Example CN Usage
cn.gbm <- pull.cn.data.TCGA(project = "TCGA-GBM", save.filename = "/home/evem/CN_Seg/Outputs")
cn.BRCA <- pull.cn.data.TCGA(project = "TCGA-BRCA", save.filename = "/home/evem/CN_Seg/Outputs")
cn.LGG <- pull.cn.data.TCGA(project = "TCGA-LGG", save.filename = "/home/evem/CN_Seg/Outputs")


# Pulling TCGA mutation data and converting to a GRanges object
pull.mutation.data.TCGA <- function(project = "", workflow.type = "STAR - Counts", save.filename = ""){
  query.mutation <- GDCquery(
    project = project,
    data.category = "Simple Nucleotide Variation",
    access = "open",
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  GDCdownload(query.mutation)
  mutation.data <- GDCprepare(query.mutation) 
  return(mutation.data)  
}
# Converting to GRanges object from DataFrame
makeGRangesFromDataFrame(mutation.gbm,
                         keep.extra.columns = TRUE,
                         ignore.strand = TRUE,
                         seqinfo = NULL,
                         seqnames.field = "Chromosome",
                         start.field = c("Start", "Start_Position"),
                         end.field = c("End", "End_Position"),
                         strand.field = "Strand") -> mutation_grange.gbm

# Example mutation usage
mutation.gbm <- pull.mutation.data.TCGA(project = "TCGA-GBM", save.filename = "/home/evem/CN_Seg/Outputs")

# Creating a GRanges object from TCGA expression data
# Function to download expression data
pull.expression.data.TCGA <- function(project = "", workflow.type = "STAR - Counts", gtf.file = ""){
  # Query and download the expression data
  query.exp <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = workflow.type)
  GDCdownload(query.exp)
  expression.data <- GDCprepare(query = query.exp)
  
  return(expression.data)
}
 
# Function to create a GRanges object from the expression data
create.expression.grange.from.TCGA <- function(expression.data, gtf.file = "") {

  # Create a TxDb object from the GTF file
  gencode.db <- makeTxDbFromGFF(gtf.file, format="gtf")
  
  # Get the gene coordinates
  gene.coords <- genes(gencode.db)
  
  # Match the gene IDs in the expression data with the gene names in the gene coordinates
  gene.ids <- rownames(assay(expression.data))
  index <- match(gene.ids, names(gene.coords))
  index <- index[!is.na(index)]
  
  # Create a GRanges object with the matched gene coordinates
  exp.grange <- gene.coords[index]
  exp.grange <- exp.grange[!is.na(names(exp.grange))]
  unique.gene.ids <- names(gene.coords)[index]
  names(exp.grange) <- unique.gene.ids
  
  return(exp.grange)
}

#### Example usage for expression data ####
expression.gbm <- pull.expression.data.TCGA(project = "TCGA-GBM", gtf.file = "/home/evem/CN_Seg/Archive/gencode.v25lift37.annotation.gtf")
exp.grange.gbm <- create.expression.grange.from.TCGA(expression.gbm, gtf.file = "/home/evem/CN_Seg/Archive/gencode.v25lift37.annotation.gtf")
expression.BRCA <- pull.expression.data.TCGA(project = "TCGA-BRCA", gtf.file = "/home/evem/CN_Seg/Archive/gencode.v25lift37.annotation.gtf")
exp.grange.BRCA <- create.expression.grange.from.TCGA(expression.BRCA, gtf.file = "/home/evem/CN_Seg/Archive/gencode.v25lift37.annotation.gtf")

# Creating a GRanges object from TCGA methylation data
# Function to download methylation data
pull.methylation.data.TCGA <- function(project = "", workflow.type = "STAR - Counts", gtf.file = ""){
  # Query and download the methylation data
  query.meth <- GDCquery(
    project = project,
    data.category = "DNA Methylation",
    data.type = "Methylation Beta Value",
    platform = "Illumina Human Methylation 450")
  GDCdownload(query.meth)
  methylation.data <- GDCprepare(query = query.meth)
  
  return(methylation.data)
}
  
# Function to create a GRanges object from the methylation data
create.meth.grange.from.TCGA <- function(methylation.data, gtf.file = "") {
  # Create a TxDb object from the GTF file
  gencode.db <- makeTxDbFromGFF(gtf.file, format="gtf")
  
  # Get the gene coordinates
  gene.coords <- genes(gencode.db)
  
  # Match the gene IDs in the methylation data with the gene names in the gene coordinates
  gene.ids <- rownames(assay(methylation.data))
  index <- match(gene.ids, names(gene.coords))
  index <- index[!is.na(index)]
  
  # Create a GRanges object with the matched gene coordinates
  meth.grange <- gene.coords[index]
  meth.grange <- meth.grange[!is.na(names(meth.grange))]
  unique.gene.ids <- names(gene.coords)[index]
  names(meth.grange) <- unique.gene.ids
  
  return(meth.grange)
}

# Example usage
methylation.gbm <- pull.methylation.data.TCGA(project = "TCGA-GBM", gtf.file = "/home/evem/CN_Seg/Archive/gencode.v25lift37.annotation.gtf")
meth.grange.gbm <- create.meth.grange.from.TCGA(methylation.gbm, gtf.file = "/home/evem/CN_Seg/Archive/gencode.v25lift37.annotation.gtf")

# Downloading clinical data
pull.clinical.data.TCGA <- function(project = "", workflow.type = "STAR - Counts", save.filename = "") {
  query.clinical <- GDCquery(
    project = project,
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    data.format = "bcr xml")
  GDCdownload(query.clinical)
  clinical.data <- GDCprepare_clinic(query = query.clinical, clinical.info = "patient")
  
  return(clinical.data)
}

# Example usage 
clinical.gbm <- pull.clinical.data.TCGA(project = "TCGA-GBM", save.filename = "/home/evem/CN_Seg/Outputs")

# Checking if column order matches across GRanges objects
check_column_order <- function(exp.grange.gbm, cn_grange.gbm) {
  # Extract column names from the metadata of the GRanges objects
  colnames1 <- colnames(mcols(exp.grange.gbm))
  colnames2 <- colnames(mcols(cn_grange.gbm))
  
  # Check if the column names are the same and in the same order
  identical(colnames1, colnames2)
}

# Example usage
exp.grange.gbm <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr2")),
                   ranges = IRanges(1:4, width = 1:4),
                   strand = strand(c("-", "+", "-", "+")),
                   score = 1:4,
                   GC = seq(1, 2, length.out = 4))

cn_grange.gbm <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr2")),
                   ranges = IRanges(1:4, width = 1:4),
                   strand = strand(c("-", "+", "-", "+")),
                   score = 1:4,
                   GC = seq(1, 2, length.out = 4))

check_column_order(exp.grange.gbm, cn_grange.gbm)  # Returns TRUE if column order matches, FALSE otherwise

