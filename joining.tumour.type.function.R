#Function for joining tumour types prior to filtering based on gene expression
#See if this effects the models in comparison to filtering separatly then merging
join.filter.tumour.types <- function(project_name1, project_name2) {
  # Construct the variable names
  exp_data_name1 <- paste0("matched.", project_name1, ".exp.assay.data.df")
  exp_data_name2 <- paste0("matched.", project_name2, ".exp.assay.data.df")
  cn_data_name1 <- paste0("matched.", project_name1, ".seg.means.df")
  cn_data_name2 <- paste0("matched.", project_name2, ".seg.means.df")
  
  # Get the data from the global environment
  exp_data1 <- get(exp_data_name1, envir = .GlobalEnv)
  exp_data2 <- get(exp_data_name2, envir = .GlobalEnv)
  cn_data1 <- get(cn_data_name1, envir = .GlobalEnv)
  cn_data2 <- get(cn_data_name2, envir = .GlobalEnv)
  
  # Combine the data
  joined.exp.data <- cbind(exp_data1, exp_data2)
  joined.cn.data <- cbind(cn_data1, cn_data2)
  joined.tumour.type <- c(rep(project_name1, ncol(exp_data1)), 
                          rep(project_name2, ncol(exp_data2)))
  
  # Filter the genes
  index.genes <- apply(joined.exp.data, 1, gp.style.filter, fold.change=3, delta=10, 
                       prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE)
  filt.joined.exp.data <- joined.exp.data[index.genes,]
  filt.joined.cn.data <- joined.cn.data[index.genes,]
  
  # Assign the results to the global environment
  assign(paste0("filt.joined.exp.data.", project_name1, ".", project_name2), 
         filt.joined.exp.data, envir = .GlobalEnv)
  assign(paste0("filt.joined.cn.data.", project_name1, ".", project_name2), 
         filt.joined.cn.data, envir = .GlobalEnv)
  assign(paste0("joined.tumour.type.", project_name1, ".", project_name2), 
         joined.tumour.type, envir = .GlobalEnv)
  

}

join.filter.tumour.types("ESCA", "ACC")

exp.ESCA.grange.2.df.matched <- exp.ESCA.grange.2.df[index.genes,]
assign(paste0("exp.ESCA.grange.2.df.matched"))
################################################################################

#Code for merging when each tumour type has seperately been filtered
#Finding common genes for when using additional tumour type as covariate
common_genes <- intersect(rownames(filt.matched.ACC.seg.means.df), 
                          rownames(filt.matched.ESCA.seg.means.df))
#Subset dataframes to only include common genes
filt.matched.ACC.seg.means.df.common <- filt.matched.ACC.seg.means.df[common_genes, ]
filt.matched.ESCA.seg.means.df.common <- filt.matched.ESCA.seg.means.df[common_genes, ]
filt.matched.ESCA.exp.assay.data.df.common <- filt.matched.ESCA.exp.assay.data.df[common_genes, ]
filt.matched.ACC.exp.assay.data.df.common <- filt.matched.ACC.exp.assay.data.df[common_genes, ]

#Checking dimensions of the two tumour types data
#Should now have the same number of rows(genes) but different number of columns(samples)
dim(filt.matched.ACC.seg.means.df.common)
dim(filt.matched.ESCA.seg.means.df.common)

#Binds 2 tumour type expression assay data together
ESCA.ACC.joined.data.exp <- cbind(filt.matched.ACC.exp.assay.data.df.common, filt.matched.ESCA.exp.assay.data.df.common)
#Binds 2 tumour type copy number data together 
ESCA.ACC.joined.data.cn <- cbind(filt.matched.ACC.seg.means.df.common, filt.matched.ESCA.seg.means.df.common)
#Creates list of values indicating what tumour type the sample is
ESCA.ACC.joined.tumour.type <- c(rep("ACC", ncol(filt.matched.ACC.exp.assay.data.df.common)), 
                                 rep("ESCA", ncol(filt.matched.ESCA.exp.assay.data.df.common)))
