#Package for GAMs
library(mgcv)
#Package for formating data in a table
library(knitr)
#Package for foreach command
library(foreach)

#Function for creating the GAM models 
#Expression as response variable
#CN as smooth term predictor variable
run.GAM.model <- function(project_name, exp, cn, i){
  #Unlist the data
  exp[i,] -> exp
  cn[i,] -> cn
  exp.data.unlist <- unlist(exp)
  seg.means.unlist <- unlist(cn)
  
  #Combining all data to call as data = combined.df
  combined.df <- data.frame( exp.data.unlist,
                             seg.means.unlist)
  #fit the model
  GAM.model <- gam(exp.data.unlist ~ s(seg.means.unlist), data = combined.df, 
                   method = "REML")
  return(GAM.model)
}

################################################################################

#Function for creating relevant outputs of GAMs
run.GAM <- function(project_name, exp, cn, i){
  #Unlist the data
  exp.data.unlist <- unlist(exp)
  seg.means.unlist <- unlist(cn)
  
  #Catch any errors
  tryCatch({
    
  #fit the model
  GAM.model <- gam(exp.data.unlist ~ s(seg.means.unlist))
  
  #Get the summary
  temp.summary <- summary(GAM.model)
  temp.p.table <- temp.summary$p.table
  temp.s.table <- temp.summary$s.table

  
  #Create output data frame
  out.df <- c(intercept.coeff = temp.summary$p.coeff,
              edf = temp.s.table[1, "edf"], reference.df = temp.s.table[1, "Ref.df"],
              F.statistic = temp.s.table[1, "F"], p.val.F = temp.s.table[1, "p-value"],
              AIC = GAM.model$aic, deviance.explained = temp.summary$dev.expl)
  
  return(out.df)
  }, error = function(e) {
  #If an error occurs it will print the error message and returns NULL
    print(paste("Error in row", project_name, ":", e$message))
    return(rep(NA,7))
})
}

#Putting rownames as ensembl gene IDs
rownames(ACC_GAM_results) <- rownames(filt.matched.ACC.exp.assay.data.df)

################################################################################

# Adjusting the p-values
smooth.term.adj.p <- p.adjust(ACC_GAM_results[, "smooth.term.p"], method = "BH")

# Add the adjusted p-values back into your results
ACC_GAM_results$smooth.term.adj.p <- smooth.term.adj.p

################################################################################

#Save the results to a file
write.csv(ACC_GAM_results, file = "/home/evem/CN_Seg/GAM_outputs/ACC_GAM_outputs.csv")
#Turn csv file into data frame then easy to interpret table
ACC_GAM_results_table <- read.csv("/home/evem/CN_Seg/GAM_outputs/ACC_GAM_outputs.csv")
ACC_GAM_results_table <- kable(ACC_GAM_results_table)
#Save GAM table to a file
sink("/home/evem/CN_Seg/GAM_outputs/ACC_GAM_results_table")
print(kable(ACC_GAM_results_table))
sink()

################################################################################

#Adding GRange data as metadata columns to GAM results
match(rownames(ACC_GAM_results), exp.ACC.grange.2$gene_id) -> ID.index
exp.ACC.grange.2.match <- exp.ACC.grange.2[ID.index]
exp.ACC.grange.2.match -> exp.ACC.grange.gam.results
cbind(mcols(exp.ACC.grange.2.match), ACC_GAM_results) -> mcols(exp.ACC.grange.gam.results)

################################################################################

data.frame(LGG_GAM_results) -> LGG_GAM_results.df

head(LGG_GAM_results.df[order(LGG_GAM_results.df$F.statistic, decreasing = T),],10)

################################################################################





GAM.model.plot <- function(gam.model, pdf.out = NULL){
  if(!is.null(pdf.out)){
  pdf(pdf.out)
  }
  ### extract the bits from the model and draw all the nice graphs you want to make
  plot(gam.model)
  ### add them here
  
  if(!is.null(pdf.out)){
    dev.off()
  }
}






################################################################################

hist(LGG_GAM_results[,"smooth.term.p"])
hist(p.adjust(LGG_GAM_results[,"smooth.term.p"], "BH"))
p.adjust(LGG_GAM_results$smooth.term.p, "BH") -> LGG_GAM_results$smooth.term.adj.p

exp.LGG.grange -> exp.LGG.grange.gam.results
cbind(mcols(exp.LGG.grange), LGG_GAM_results) -> mcols(exp.LGG.grange.gam.results)

data.frame(LGG_GAM_results) -> LGG_GAM_results.df

head(LGG_GAM_results.df[order(LGG_GAM_results.df$F.statistic, decreasing = T),],10)




#Save the results to a file
write.csv(LGG_GAM_results, file = "/home/evem/CN_Seg/GAM_outputs/LGG_GAM_outputs.csv")

#Turn csv file into data frame then easy to interpret table
LGG_GAM_results_table <- read.csv("/home/evem/CN_Seg/GAM_outputs/LGG_GAM_outputs.csv")
LGG_GAM_results_table <- kable(LGG_GAM_results_table)
#Save GAM table to a file
sink("/home/evem/CN_Seg/GAM_outputs/LGG_GAM_results_table")
print(kable(LGG_GAM_results_table))
sink()


#GBM data
#Open a PDF file to save the residual plots
#pdf(file = paste0("/home/evem/CN_Seg/GAM_outputs", "/GAM_Residuals_GBM.pdf"))

#Check number of rows for calling GAM function
nrow(filt.matched.GBM.exp.assay.data.df)

#Call run.GAM function row-wise (per gene)
GBM_GAM_results <- foreach(i = 1:10105, .combine = rbind)%do%{
  run.GAM(project_name = "GBM", filt.matched.GBM.exp.assay.data.df[i,],filt.matched.GBM.seg.means.df[i,], i)}

#Close the PDF file
#dev.off()

#Save the results to a file
write.csv(GBM_GAM_results, file = "/home/evem/CN_Seg/GAM_outputs/GBM_GAM_outputs.csv")

#Turn csv file into data frame then easy to interpret table
GBM_GAM_results_table <- read.csv("/home/evem/CN_Seg/GAM_outputs/GBM_GAM_outputs.csv")
GBM_GAM_results_table <- kable(GBM_GAM_results_table)
#Save GAM table to a file
sink("/home/evem/CN_Seg/GAM_outputs/GBM_GAM_results_table")
print(kable(GBM_GAM_results_table))
sink()
#############################################################################################################
