#Running linear model with an interaction
#Uses GAM argument
#Expression is predictor variable
#CN is response variable
#Covariate is interaction term along with CN
run.LM.interaction.model <- function(project_name, exp, cn, covariate, i){
  
  #Unlist the data
  exp[i,] -> exp
  cn[i,] -> cn
  exp.data.unlist <- unlist(exp)
  seg.means.unlist <- unlist(cn)
  covariate.unlist <- as.factor(unlist(covariate))
  
  #Combing all data to call as data = combined.df
  combined.df <- data.frame( exp.data.unlist,
                             seg.means.unlist,
                             covariate.unlist)
  #fit the model
  #this formula models the relationship between CN and expression for each level
  #of the covariate.unlist (each tumour type)
  LM.interaction.model <- lm(exp.data.unlist ~ seg.means.unlist * covariate.unlist,
                             data = combined.df)
  return(LM.interaction.model)
}

run.LM.interaction.model(project_name = "ESCA", filt.joined.exp.data.ESCA.ACC[2,],
                          filt.joined.cn.data.ESCA.ACC[2,], 
                          covariate = joined.tumour.type.ESCA.ACC, 1) -> LM.interaction.model

################################################################################

#Function for creating table of LM with interaction outputs for each gene 
run.LM.interaction <- function(project_name, exp, cn, covariate, i){
  
  
  #Unlist the data
  exp.data.unlist <- unlist(exp)
  seg.means.unlist <- unlist(cn)
  covariate.unlist <- as.factor(covariate)
  
  #Combing all data to call as data = combined.df
  combined.df <- data.frame( exp.data.unlist,
                             seg.means.unlist,
                             covariate.unlist)
  
  #Catch any errors
  tryCatch({
    
    library(dplyr)
    
    #fit the model
    LM.interaction.model <- lm(exp.data.unlist ~ seg.means.unlist * covariate.unlist,
                                data = combined.df)
    
    #Get the summary
    temp.summary <- summary(LM.interaction.model)
    temp.coefficients <- temp.summary$coefficients
    temp.r.squared <- temp.summary$r.squared
    temp.AIC <- AIC(LM.interaction.model)
    
    
    #Create output data frame
    out.df <- c(Intercept.coef = temp.coefficients["(Intercept)", "Estimate"], 
                seg.means.reflevel.coef = temp.coefficients["seg.means.unlist", "Estimate"],
                interaction.coef = temp.coefficients["seg.means.unlist:covariate.unlistESCA", "Estimate"], 
                reflevel.coef.Std.Error = temp.coefficients["seg.means.unlist", "Std. Error"],
                interaction.coef.Std.Error = temp.coefficients["seg.means.unlist:covariate.unlistESCA", "Std. Error"],
                reflevel.coef.t.value = temp.coefficients["seg.means.unlist", "t value"],
                interaction.coef.t.value = temp.coefficients["seg.means.unlist:covariate.unlistESCA", "t value"],
                reflevel.coef.p.value = temp.coefficients["seg.means.unlist", "Pr(>|t|)"],
                interaction.coef.p.value = temp.coefficients["seg.means.unlist:covariate.unlistESCA", "Pr(>|t|)"],
                AIC = temp.AIC, r.squared = temp.r.squared)
    
    
    return(out.df)
  }, error = function(e) {
    #If an error occurs it will print the error message and returns NULL
    print(paste("Error in row", project_name, ":", e$message))
    return(rep(NA,11))
  })
}

#Call run.LM.interaction function row-wise (per gene)
#Creates output with desired model metrics in table
ESCA_interaction_LM_results <- foreach(i = 1:nrow(filt.joined.exp.data.ESCA.ACC), .combine = rbind)%do%{
  run.LM.interaction(project_name = "ESCA", exp = filt.joined.exp.data.ESCA.ACC[i,], 
                      cn = filt.joined.cn.data.ESCA.ACC[i,], 
                      covariate = joined.tumour.type.ESCA.ACC)}

#Putting rownames as gene names
rownames(ESCA_interaction_LM_results) <- rownames(filt.joined.exp.data.ESCA.ACC)
