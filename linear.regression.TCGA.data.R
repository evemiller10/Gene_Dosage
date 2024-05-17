#Package for GAMs
library(mgcv)
#Package for formating data in a table
library(knitr)
#Package for foreach command
library(foreach)

#Using GAM argument to create linear model
#Function for creating the GAM models 
run.LM.model <- function(project_name, exp, cn, i){
  #Unlist the data
  exp[i,] -> exp
  cn[i,] -> cn
  exp.data.unlist <- unlist(exp)
  seg.means.unlist <- unlist(cn)
  #fit the model
  LM.model <- gam(exp.data.unlist ~ seg.means.unlist)
  return(LM.model)
}

run.LM.model(project_name = "ACC", exp = filt.matched.ACC.exp.assay.data.df["ENSG00000001036.14", ],
             cn = filt.matched.ACC.seg.means.df["ENSG00000001036.14", ], 1) -> LM.test.model

plot(LM.test.model)


###############################################################################
run.LM <- function(project_name, exp, cn, i){
  #Unlist the data
  exp.data.unlist <- unlist(exp)
  seg.means.unlist <- unlist(cn)
  
  #Catch any errors
  tryCatch({
    
    #fit the model
    LM.model <- gam(exp.data.unlist ~ seg.means.unlist)
    
    #Get the summary
    temp.summary <- summary(LM.model)
    temp.coefficients <- temp.summary$coefficients
    
    
    #Create output data frame
    out.df <- c(Int.coef = temp.summary$p.coeff["(Intercept)"], 
                coef = temp.summary$p.coeff["seg.means.unlist"])
    
    
    return(out.df)
  }, error = function(e) {
    #If an error occurs it will print the error message and returns NULL
    print(paste("Error in row", project_name, ":", e$message))
    return(rep(NA,2))
  })
}

nrow(filt.matched.ACC.exp.assay.data.df)
#Call run.GAM function row-wise (per gene)
ACC_LM_results <- foreach(i = 1:nrow(filt.matched.ACC.exp.assay.data.df), .combine = rbind)%do%{
  run.LM(project_name = "ACC", filt.matched.ACC.exp.assay.data.df[i,],
         filt.matched.ACC.seg.means.df[i,], i)}

rownames(ACC_LM_results) <- rownames(filt.matched.ACC.exp.assay.data.df)

data.frame(ACC_LM_results) -> ACC_LM_results.df

head(ACC_LM_results.df[order(ACC_LM_results.df$coef.seg.means.unlist, decreasing = F),],10)

#################################################################################

#Adding adjusted p-values to results
# Adjusting the p-values
adj.p.value.intercept <- p.adjust(ACC_GAM_results[, "p.val.intercept"], method = "BH")
adj.p.value.slope <- p.adjust(ACC_GAM_results[, "p.val.slope"], method = "BH")

# Add the adjusted p-values back into your results
ACC_LM_results$adj.p.value.intercept <- adj.p.value.intercept
ACC_LM_results$adj.p.value.slope <- adj.p.value.slope

################################################################################

#Save the results to a file
write.csv(ACC_LM_results, file = "/home/evem/CN_Seg/LM_outputs/ACC_LM_outputs.csv")
#Turn csv file into data frame then easy to interpret table
ACC_LM_results_table <- read.csv("/home/evem/CN_Seg/LM_outputs/ACC_LM_outputs.csv")
ACC_LM_results_table <- kable(ACC_LM_results_table)
#Save GAM table to a file
sink("/home/evem/CN_Seg/LM_outputs/ACC_LM_results_table")
print(kable(ACC_LM_results_table))
sink()

################################################################################
#Function for creating the Linear models 
run.LM.model <- function(project_name, x, y, i){
  #Unlist the data
  x[i,] -> x
  y[i,] -> y
  exp.data.unlist <- unlist(x)
  seg.means.unlist <- unlist(y)
  #fit the model
  LM.model <- lm(exp.data.unlist ~ seg.means.unlist)
  return(LM.model)
  
}

run.LM.model(project_name = "ACC", x = filt.matched.ACC.exp.assay.data.df["ENSG00000001036.14",],
             y = filt.matched.ACC.seg.means.df["ENSG00000001036.14",], 1) -> LM.test.model.x

run.LM <- function(project_name, x, y, i){
  #Unlist the data
  exp.data.unlist <- unlist(x)
  seg.means.unlist <- unlist(y)
  
  #Catch any errors
  tryCatch({
    
    #fit the model
    LM.model <- lm(exp.data.unlist ~ seg.means.unlist)
    
    #Get the summary
    temp.summary <- summary(GAM.model)
    temp.p.table <- temp.summary$p.table
    temp.s.table <- temp.summary$s.table
    
    
    #Create output data frame
    out.df <- c(seg.means.coeff = GAM.model$coefficients[,2], 
                intercept.coeff = GAM.model$coefficients[,1])
    
    return(out.df)
  }, error = function(e) {
    #If an error occurs it will print the error message and returns NULL
    print(paste("Error in row", project_name, ":", e$message))
    return(rep(NA,2))
  })
}

ACC_LM_results <- foreach(i = 1:nrow(filt.matched.ACC.exp.assay.data.df), .combine = rbind)%do%{
  run.LM(project_name = "ACC", filt.matched.ACC.exp.assay.data.df[i,],
         filt.matched.ACC.seg.means.df[i,], i)}
