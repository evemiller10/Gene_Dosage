library(TCGAbiolinks)
library(SummarizedExperiment)
library(GenomicRanges)
library(GenomicFeatures)


# Function to pull in TCGA data
pull.data.TCGA <- function(
    project = "",
    workflow.type = "STAR - Counts",
    save.filename = ""
){
  
  # project = "TCGA-GBM"
  # save.filename = "/home/evem/CN_Seg/Outputs"
  # workflow.type = "STAR - Counts"
  
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
  
  # Downloading Methylation Data
  query.meth <- GDCquery(
    project = project,
    data.category = "DNA Methylation",
    data.type = "Methylation Beta Value",
    platform = "Illumina Human Methylation 450")
  GDCdownload(query.meth)
  methylation <- GDCprepare(query.meth)
  
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
  
  # Downloading Proteomics Data
  query.proteomics <- GDCquery(
    project = project,
    data.category = "Proteome Profiling",
    data.type = "Protein Expression Quantification")
  GDCdownload(query.proteomics)
  proteomics <- GDCprepare(query = query.proteomics)
  
  
  # Returns the 6 elements as a list
  return(list(expression = expression,
              clinical = clinical,
              methylation = methylation,
              copynumber = copynumber,
              mutational = mutational,
              proteomics = proteomics))
}

TCGA.GBMdata <- pull.data.TCGA(project = "TCGA-GBM", save.filename = "/home/evem/CN_Seg/Outputs")

# error message here 
TCGA.GBMdata$expression



# Converting to GRanges object from DataFrame  
convert_dataframe_to_GRanges <- function(df = "", seqnames.field = "Chromosome", start.field = c("Start", "Start_Position"), end.field = c("End", "End_Position"), strand.field = "Strand") {
  granges <- makeGRangesFromDataFrame(df,
                           keep.extra.columns = TRUE,
                           ignore.strand = TRUE,
                           seqinfo = NULL,
                           seqnames.field = seqnames.field,
                           start.field = start.field,
                           end.field = end.field,
                           strand.field = strand.field)

}
cn.gbm.grange <- convert_dataframe_to_GRanges(df = TCGA.GBMdata$copynumber)

# Function to create a GRanges object from the expression data
create.expression.grange.from.TCGA <- function(expression_data, gtf.file = "") {
  
  # Create a TxDb object from the GTF file
  gencode.db <- makeTxDbFromGFF(gtf.file, format="gtf")
  
  # Get the gene coordinates
  gene.coords <- genes(gencode.db)
  
  # Match the gene IDs in the expression data with the gene names in the gene coordinates
  gene.ids <- rownames(assay(expression_data))
  index <- match(gene.ids, names(gene.coords))
  index <- index[!is.na(index)]
  
  # Create a GRanges object with the matched gene coordinates
  exp.grange <- gene.coords[index]
  exp.grange <- exp.grange[!is.na(names(exp.grange))]
  unique.gene.ids <- names(gene.coords)[index]
  names(exp.grange) <- unique.gene.ids
  
  return(exp.grange)
}

assays(TCGA.GBMdata$expression)$"tpm_unstrand" -> expression_data_gbm
exp.gene.ids <- row_data$gene_id

# Replace 'gene_id_column' with the name of the column that contains the gene IDs
rownames(TCGA.GBMdata$expression) <- TCGA.GBMdata$expression$gene_id_column
exp.grange.gbm <- create.expression.grange.from.TCGA(expression_data_gbm, gtf.file = "/home/evem/CN_Seg/Archive/gencode.v25lift37.annotation.gtf")


# Remove metadata columns
mcols(cn.gbm.grange)$Copy_Number <- NULL
mcols(cn.gbm.grange)$Major_Copy_Number <- NULL
mcols(cn.gbm.grange)$Minor_Copy_Number <- NULL

# Create a vector of 'Call' values based on 'Seg.CN' values
call_values <- ifelse(cn.gbm.grange$Seg.CN > 0.2, "gain",
                      ifelse(cn.gbm.grange$Seg.CN < -0.2, "loss", "Normal"))
# Add 'Call' metadata column with specific values
mcols(cn.gbm.grange)$CNcall <- call_values

exp.grange.gbm.list <- GRangesList(exp.grange.gbm)

# creates a granges object with chromosomes binned at the sizes of 1Mb
chromosome.map <- read.delim("./Archive/chromosome.map")

#creates a granges object with certain custom regions of interest
#keep.extra.columns = TRUE tells the function to keep all columns from the data frame that aren't used to construct the data frame
#ignore.strand= TRUE tells the function to ignore the strand information when creating granges object
select.features <- read.delim("./Archive/select.features.txt")
makeGRangesFromDataFrame(chromosome.map, keep.extra.columns = TRUE,ignore.strand=TRUE) -> chr.grange
makeGRangesFromDataFrame(select.features, keep.extra.columns = TRUE,ignore.strand=TRUE) -> select.grange

#generates a list of matrices where individuals represented by columns and genes as rows
#contains 4 slots:
#seg.means= average CN ration within each segment
#CN
#CNcall= categorical calls representing aplifications/deletions etc
#num.mark= number of markers within each segment
segMatrix(exp.grange.gbm, cn.gbm.grange) -> seg.list.gbm
segMatrix(chr.grange,cn.gbm.grange) -> seg.chr.list
segMatrix(select.grange,cn.gbm.grange) -> seg.select.list

#creates colour palette 
#n=100 sets the split of colours in the gradient
colors <- colorRampPalette(c('darkblue','blue','white','red','darkred'))(n=100)
library(NMF)
# Remove NA values from each element of the list
seg.select.list$seg.means <- na.omit(seg.select.list$seg.means)
seg.select.list$CNcall <- na.omit(seg.select.list$CNcall)
seg.select.list$num.mark <- na.omit(seg.select.list$num.mark)

gsub("\\..*","",rownames(seg.list.gbm$CNcall)) -> rownames(seg.list.gbm$CNcall)
gsub("\\..*","",rownames(TCGA.GBMdata$expression)) -> rownames(TCGA.GBMdata$expression)
gsub("\\..*","",names(exp.grange.gbm)) -> names(exp.grange.gbm)



#match the expression data column IDs to the available copy number sample IDs
match(colnames(seg.list.gbm$CNcall), colnames(TCGA.GBMdata$expression)) -> index
TCGA.GBMdata$expression[,index[!is.na(index)]] -> select.exp.data

#match the copy number sample IDs to the available expression data column IDs
#this is where error was occuring as seg.list$CNcall was previously not a matrix
seg.list.gbm$CNcall[,which(!is.na(index))] -> matched.cn
seg.list.gbm$seg.means[,which(!is.na(index))] -> matched.seg.means

min.number.samples <- 5


# set up normal margins
# the calculate normal
normal.values <- calculate.normal.vals(matched.cn, select.exp.data, cores =20)
sample.number <- 10






