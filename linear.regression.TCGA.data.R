#Package for GAMs
library(mgcv)
#Package for formating data in a table
library(knitr)
#Package for foreach command
library(foreach)

#Using GAM argument to create linear model
#Function for running a linear model
run.LM.model <- function(project_name, exp, cn, i){
  
  #Unlist the data
  exp[i,] -> exp
  cn[i,] -> cn
  exp.data.unlist <- unlist(exp)
  seg.means.unlist <- unlist(cn)
  
  #Combining data so can return a data frame of expression and CN data
  #Needed for plotting of LM
  combined.df <- data.frame( exp.data.unlist,
                             seg.means.unlist)
  
  #fit the model
  LM.model <- gam(exp.data.unlist ~ seg.means.unlist)
  return(list(model = LM.model, data = data.frame(seg.means.unlist,
                                                  exp.data.unlist)))
}

#Call for running ACC LM model
run.LM.model(project_name = "ACC", exp = filt.matched.ACC.exp.assay.data.df[2, ],
             cn = filt.matched.ACC.seg.means.df[2, ], 1) -> LM.test.model

#Plotting the linear model
#No automatic plot function as for gam with smooth functions
#Have to return expression and CN data as dataframes then plot using these
plot(LM.test.model$data$seg.means.unlist, LM.test.model$data$exp.data.unlist, 
     xlab = "copy number", ylab = "expression")

# Add the regression line to the plot
abline(LM.test.model$model, col = "red")

##############################################################################################

#Function for creating table of LM outputs for each gene in selected tumour type
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
    out.df <- c(Integer.coef = temp.summary$p.table["(Intercept)", "Estimate"], 
                seg.means.coef = temp.summary$p.table["seg.means.unlist", "Estimate"], 
                coef.Std.Error = temp.summary$p.table["seg.means.unlist", "Std. Error"],
                coef.t.value = temp.summary$p.table["seg.means.unlist", "t value"],
                coef.p.value = temp.summary$p.table["seg.means.unlist", "Pr(>|t|)"],
                AIC = LM.model$aic, deviance.explained = temp.summary$dev.expl)
    
    
    return(out.df)
  }, error = function(e) {
    #If an error occurs it will print the error message and returns NULL
    print(paste("Error in row", project_name, ":", e$message))
    return(rep(NA,7))
  })
}

###########################################################################################

#Adding adjusted p-values to results
# Adjusting the p-values
adj.p.value.intercept <- p.adjust(ACC_GAM_results[, "p.val.intercept"], method = "BH")
adj.p.value.slope <- p.adjust(ACC_GAM_results[, "p.val.slope"], method = "BH")

# Add the adjusted p-values back into your results
ACC_LM_results$adj.p.value.intercept <- adj.p.value.intercept
ACC_LM_results$adj.p.value.slope <- adj.p.value.slope

###########################################################################################

#Save the results to a file
write.csv(ACC_LM_results, file = "/home/evem/CN_Seg/LM_outputs/ACC_LM_outputs.csv")
#Turn csv file into data frame then easy to interpret table
ACC_LM_results_table <- read.csv("/home/evem/CN_Seg/LM_outputs/ACC_LM_outputs.csv")
ACC_LM_results_table <- kable(ACC_LM_results_table)
#Save GAM table to a file
sink("/home/evem/CN_Seg/LM_outputs/ACC_LM_results_table")
print(kable(ACC_LM_results_table))
sink()
