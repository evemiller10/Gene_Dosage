source ('./Scripts/exp2cn.source.R')


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicFeatures")
library(GenomicRanges)
library(GenomicFeatures)
install.packages("AnnotationDbi")
library(AnnotationDbi)
library(doMC)

# converts segmentation file created in crllm/cbs into a genomic ranges object
# the seg file is read into a genomicranges object 
# the chromosomes are given in seqnames, the start and end as the ranges input
# individual samples are given in the extra columns under header sample
# Num.markers indicates the number of data points within the segment
# Seg.CN indicates log ratio where 0 is normal <0 indicates potential loss and >0 indicates potential gain and call is the significance call associated with that segment
readSegFile("./Archive/weighted.smoothed.seg") -> seg.grange

# Renaming Call metadata column as CNcall
names(mcols(seg.grange))[names(mcols(seg.grange)) == 'Call'] <- 'CNcall'

# gencode.db <- makeTxDbFromGFF("/home/evem/gencode.v25lift37.annotation.gtf", format="gtf")

# Save the transcriptDb object as a sql database object
# saveDb (gencode.db, file = "/home/evem/gencode_human_v25.sqlite")

# creates an expression genomic ranges list from the gencode data used to create the RNA data
# each slot in the list represents an individual ensembl gene ENSG is an ensembl ID
exp.grange <- readGencodeDb("./Archive/gencode_human_v25.sqlite")

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
gencode.db <- makeTxDbFromGFF("/home/evem/CN_Seg/Archive/gencode.v25lift37.annotation.gtf", format="gtf")
segMatrix(exp.grange,seg.grange) -> seg.list
segMatrix(chr.grange,seg.grange) -> seg.chr.list
segMatrix(select.grange,seg.grange) -> seg.select.list

#save(seg.chr.list, file="/home/evem/CN_Seg/Outputs/seg.chr.list.Rdata")
#save(seg.select.list, file="/home/evem/CN_Seg/Outputs/seg.select.list.Rdata")
#save(seg.list, file="/home/evem/CN_Seg/Outputs/seg.list.Rdata")
#creates colour palette 
#n=100 sets the split of colours in the gradient
colors <- colorRampPalette(c('darkblue','blue','white','red','darkred'))(n=100)
library(NMF)
aheatmap(as.matrix(seg.select.list[[1]]), Rowv=NA, color=colors, breaks = 0)

###############################################################################

# reads in expression data matrix, rownames are genes and column names are sample IDs
exp.data <- read.delim(file = "./Archive/vsd.gds.primary.txt")

#remove ensembl id .ext
#replaces with empty string
#turning seg.list$CNcall from characeter vector to matrix
gsub("\\..*","",rownames(seg.list$CNcall)) -> rownames(seg.list$CNcall)
gsub("\\..*","",rownames(exp.data)) -> rownames(exp.data)
gsub("\\..*","",names(exp.grange)) -> names(exp.grange)

#create an index select and match genes
#matches the row names of seg.list$CNcall and exp.data
#subsets data to only keep rows where the corresponding row in index is not missing
match (rownames(seg.list$CNcall), rownames(exp.data)) -> index
exp.grange[which(!is.na(index))] -> exp.grange
exp.data[index[!is.na(index)],] -> exp.data
seg.list$CNcall[which(!is.na(index)),] -> seg.list$CNcall
seg.list$seg.means[which(!is.na(index)),] -> seg.list$seg.means

#match the expression data column IDs to the available copy number sample IDs
match(colnames(seg.list$CNcall), colnames(exp.data)) -> index
exp.data[,index[!is.na(index)]] -> select.exp.data

#match the copy number sample IDs to the available expression data column IDs
#this is where error was occuring as seg.list$CNcall was previously not a matrix
seg.list$CNcall[,which(!is.na(index))] -> matched.cn
seg.list$seg.means[,which(!is.na(index))] -> matched.seg.means


# calculate a very basic relationship between copy number and expression
# set the minumum number of samples possessing copy number changes 
# before consideration for stats test
min.number.samples <- 5

# set up normal margins
# the calculate normal
normal.values <- calculate.normal.vals(matched.cn, select.exp.data, cores =20)
sample.number <- 10

# plot lowess smoothing expression and copy number chromosome by chromosome
chromosomal.plot(normal.values, 
select.exp.data, 
exp.grange, 
sample.number = sample.number, 
matched.seg.means, 
matched.cn = matched.cn, 
output.file = "Outputs/chromosome.pdf")

###########################################################################
t.test.cna(matched.cn, select.exp.data, min.number.samples) -> all.tests

# produces list detailing t statistic (differences in expression between samples with relevant copy number change vs normal sample)
# the 5 slots in the list pertain to the following types of copy number change in order:
# "focal.gain"
# "focal.amp"
# "broad.gain"
# "broad.del"
# "focal.del"
names(exp.grange) -> gene.names
as.character(seqnames(exp.grange)) -> chromosome.locs
as.numeric(start(exp.grange)) -> start.locs
levels(as.factor(chromosome.locs)) -> all.chrom

locs <- cbind(chromosome.locs,start.locs)
rownames(locs) <- gene.names

pdf("./Outputs/dosage.dependency.chr17.pdf", width=20,height=10)
all.tests[[3]] -> focal.gain.test
cbind(locs[rownames(focal.gain.test),],focal.gain.test) -> focal.gain.test
focal.gain.test[focal.gain.test[,1]=="chr17",] ->select.focal.gain.test
color <- ifelse(as.numeric(select.focal.gain.test[,5])<0.05, "orange", "yellow")
size <- ifelse(as.numeric(select.focal.gain.test[,5])<0.05,1,0.5)
plot(as.numeric(select.focal.gain.test[,2]),as.numeric(select.focal.gain.test[,3]), pch = 19, xlim = c(1,1.3E08), ylim =c(-6,21), col = color, cex = size, xlab = "bp", ylab = "CN dosage dependency")
abline(h=0)
select.focal.gain.test[color=="orange",]

all.tests[[2]] -> focal.amp.test
cbind(locs[rownames(focal.amp.test),],focal.amp.test) -> focal.amp.test
focal.amp.test[focal.amp.test[,1]=="chr17",] ->select.focal.amp.test
color <- ifelse(as.numeric(select.focal.amp.test[,5])<0.05, "darkred", "red")
size <- ifelse(as.numeric(select.focal.amp.test[,5])<0.05,1,0.5)
points(as.numeric(select.focal.amp.test[,2]),as.numeric(select.focal.amp.test[,3]), pch = 19, col = color, cex = size)
select.focal.amp.test[color=="darkred",]
legend(1,20,c("focal gain - Significant","focal gain - Non-significant","broad gain - Significant","broad gain - Non-Significant"), pch = 19, col = c("darkred","red","orange","yellow"))
dev.off()
```