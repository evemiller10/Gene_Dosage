# Singular function to pull in TCGA data
# Expression
# Copy Number
# Mutation
# Clinical
# Returns each data type as element in a list
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
  query.maf <- GDCquery(
    project = project,
    data.category = "Simple Nucleotide Variation",
    access = "open",
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  GDCdownload(query.maf)
  mutational <- GDCprepare(query.maf)
  
  # Downloading Clinical Data
  query.clin <- GDCquery(
    project = project,
    #file.type = "xml",
    data.category = "Clinical",
    #legacy = FALSE,
    data.type = "Clinical Supplement",
    data.format = "BCR XML")
  GDCdownload(query.clin)
  clinical <- GDCprepare(query = query.clin)
  
  # Returns the elements as a list
  return(list(expression = expression,
              clinical = clinical,
              copynumber = copynumber,
              mutational = mutational))
  
}


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

# Function to convert RangedSummarizedExperiment into GRanges object
rse_to_GRanges <- function(rse) {
  grange <- rowRanges(rse)
  # Add "sample' column as metadata
  mcols(grange)$sample <- colData(rse)$sample
  
  return(grange)
}

# Function for filtering genes
gp.style.filter <- function (x,fold.change, delta, prop, base, prop.base, na.rm = TRUE, neg.rm = TRUE)
{
  if (na.rm == TRUE){x <- x[!is.na(x)]}
  if (neg.rm == TRUE){x <- x[x > 0]}
  lenx <- length(x)
  if (lenx < 4 || lenx < prop + 1){return(FALSE)}
  srtd <- sort(x)
  if (prop < 1) {
    bot <- lenx * prop
  }else {bot <- prop}
  top <- lenx - bot
  fold.filter<-FALSE
  delta.filter<-FALSE
  base.filter<-FALSE
  if (max(srtd[bot:top]) / min(srtd[bot:top]) > fold.change){fold.filter<-TRUE}
  if (max(srtd[bot:top]) - min(srtd[bot:top]) > delta){delta.filter<-TRUE}
  if (length(which(srtd > base))>(prop.base*lenx)){base.filter<-TRUE}
  if (fold.filter ==TRUE & delta.filter == TRUE & base.filter == TRUE){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
