#Package for GAMs
library(mgcv)
#Package for formating data in a table
library(knitr)
#Package for foreach command
library(foreach)

#Function for running GAM model with categorical interaction term
#Built of use for interactions term being 2 tumour types
run.GAM.model.interaction <- function(project_name, exp, cn, covariate, covariate_type, i){

  # Unlist the data
  exp[i,] -> exp
  cn[i,] -> cn
  #covariate[i,] -> covariate
  exp.data.unlist <- unlist(exp)
  seg.means.unlist <- unlist(cn)
  covariate.unlist <- as.factor(unlist(covariate))
  
  #Combing all data to call as data = combined.df
  combined.df <- data.frame( exp.data.unlist,
                    seg.means.unlist,
                    covariate.unlist
                    )
  library(dplyr)

  # Fit the model
  if (covariate_type == "linear") {
    GAM.model <- gam(exp.data.unlist ~ covariate.unlist + 
                       s(seg.means.unlist, by = covariate.unlist), 
                     data = combined.df, method = "REML") 
  } else if (covariate_type == "smooth") {
    GAM.model <- gam(exp.data.unlist ~ s(seg.means.unlist, covariate.unlist), 
                     data = combined.df, method = "REML")
  } else {
    stop("Invalid covariate type. Choose 'linear' or 'smooth'.")
  }
  
  return(GAM.model)
}

#Calling the GAM with covariate model function
run.GAM.model.interaction(project_name = "ESCA", filt.joined.exp.data.ESCA.ACC[1,],
                                                    filt.joined.cn.data.ESCA.ACC[1,], 
                                         covariate = joined.tumour.type.ESCA.ACC, 
                                         covariate_type = "linear", 1) -> row.test.interaction.model

################################################################################
#Function for creating relevant outputs of GAMs with categorical intercation term
run.GAM.interaction <- function(project_name, exp, cn, covariate, covariate_type, i){
 
   #Unlist the data
  exp.data.unlist <- unlist(exp)
  seg.means.unlist <- unlist(cn)
  covariate.unlist <- as.factor(covariate)
  
  #Combing all data to call as data = combined.df
  combined.df <- data.frame( exp.data.unlist,
                             seg.means.unlist,
                             covariate.unlist
  )
  
  #Catch any errors
  tryCatch({
    
    library(dplyr)
    
    # Fit the model
    if (covariate_type == "linear") {
      GAM.model <- gam(exp.data.unlist ~ covariate.unlist + s(seg.means.unlist, 
                      by = covariate.unlist), data = combined.df, method = "REML") 
    } else if (covariate_type == "smooth") {
      GAM.model <- gam(exp.data.unlist ~ s(seg.means.unlist, covariate.unlist), 
                       data = combined.df, method = "REML")
    } else {
      stop("Invalid covariate type. Choose 'linear' or 'smooth'.")
    }
    
    
    #Get the summary
    temp.summary <- summary(GAM.model)
    temp.p.table <- temp.summary$p.table
    temp.s.table <- temp.summary$s.table
    
    #Create output data frame
    out.df <- c(GCV = GAM.model$gcv.ubre, 
                n = temp.summary$n, p.coeff = temp.p.table[, 1], std.error = temp.p.table[, 2], 
                t.val = temp.p.table[, 3], t.associated.p.value = temp.p.table[, 4],
                chi.sq = temp.summary$chi.sq,
                deviance = temp.summary$dev.expl,
                edf.1 = temp.s.table[1,1], F.statistic.1 = temp.s.table[1,3],
                smooth.term.p.1 =temp.s.table[1,4], edf.2 = temp.s.table[2,1],
               F.statistic.2 = temp.s.table[2,3], smooth.term.p.2 = temp.s.table[2,4])
    
    #Creates variables with project abbreviations in names
    #assign(paste0("GAM.model.", project_name, ".", i), GAM.model, envir = .GlobalEnv)
    
    return(out.df)
  }, error = function(e) {
    #If an error occurs it will print the error message and returns NULL
    print(paste("Error in row", project_name, ":", e$message))
    return(rep(NA,14))
  })
}


#Call run.GAM function row-wise (per gene)
ESCA_interaction_GAM_results <- foreach(i = 1:nrow(filt.joined.exp.data.ESCA.ACC), .combine = rbind)%do%{
   run.GAM.interaction(project_name = "ESCA", exp = filt.joined.exp.data.ESCA.ACC[i,], 
                       cn = filt.joined.cn.data.ESCA.ACC[i,], 
                      covariate = joined.tumour.type.ESCA.ACC, covariate_type = "linear")}

#Putting rownames as gene names
rownames(ESCA_interaction_GAM_results) <- rownames(filt.joined.exp.data.ESCA.ACC)

################################################################################

#Adding adjusted p-values to results
# Adjusting p-values
t.associated.adj.p.value <- p.adjust(ESCA_interaction_GAM_results[, "t.associated.p.value"], method = "BH")
smooth.term.adj.p.1 <- p.adjust(ESCA_interaction_GAM_results[, "smooth.term.p.1"], method = "BH")
smooth.term.adj.p.2 <- p.adjust(ESCA_interaction_GAM_results[, "smooth.term.p.2"], method = "BH")

# Add the adjusted p-values back into your results
ESCA_interaction_GAM_results$t.associated.adj.p.value <- t.associated.adj.p.value
ESCA_interaction_GAM_results$smooth.term.adj.p.1 <- smooth.term.adj.p.1
ESCA_interaction_GAM_results$smooth.term.adj.p.2 <- smooth.term.adj.p.2

###############################################################################
