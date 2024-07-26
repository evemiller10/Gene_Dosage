#Scaling gene expression values via tumour type

filt.matched.ACC.exp.assay.data.df.scale <- filt.matched.ACC.exp.assay.data.df %>%
  mutate(Tumour_Type = "ACC")


# Extract expression values (numeric matrix)
ACC_expression_matrix <- as.matrix(filt.matched.ACC.exp.assay.data.df.scale[, 
                         -ncol(filt.matched.ACC.exp.assay.data.df.scale)])

# Extract the tumour type factor
tumour_type <- filt.matched.ACC.exp.assay.data.df.scale$Tumour_Type


# Scale the expression values by tumour type
scale_centre <- function(x, tumour_type){
df <- data.frame(tumour_type = tumour_type, exp = as.numeric(x))
scaled_expression <- scale_by(exp~tumour_type, data = df)
return(as.numeric(scaled_expression))}

all_tumours_scaled_expression <- t(apply(all_tumours_joined_data_exp, 1, scale_centre, tumour_type = all_tumours_joined_tumour_type))
colnames(all_tumours_scaled_expression) <- colnames(all_tumours_joined_data_exp)

all_tumours_scaled_expression.2 <- t(apply(all.tumours.joined.data.exp, 1, scale_centre, tumour_type = all.tumours.joined.tumour.type))
colnames(all_tumours_scaled_expression.2) <- colnames(all.tumours.joined.data.exp)

# remember to try log(tpm+1) does this make a difference
# 
# scaled_expression <- foreach(i = 1:nrow(all_tumours_joined_data_exp), .combine = rbind)%do%{
#   df <- data.frame(tumour_type = all_tumours_joined_tumour_type, exp = as.numeric(all_tumours_joined_data_exp[i,]))
#   scaled_expression <- scale_by(exp~tumour_type, data = df)
#   return(as.numeric(scaled_expression))} 


boxplot(as.numeric(all_tumours_joined_data_exp[1,])~all_tumours_joined_tumour_type)
boxplot(scaled_expression[1,]~all_tumours_joined_tumour_type)

# Customizing the first boxplot (raw expression values)
boxplot(as.numeric(all_tumours_joined_data_exp[1,]) ~ all_tumours_joined_tumour_type,
        col = rainbow(length(unique(all_tumours_joined_tumour_type))), # Different colors for each tumor type
        xlab = "Tumor Type", # Custom x-axis label
        ylab = "Raw Expression Values", # Custom y-axis label
        main = "Boxplot of Raw Expression Values by Tumor Type") # Main title

# Customizing the second boxplot (scaled expression values)
boxplot(scaled_expression[1,] ~ all_tumours_joined_tumour_type,
        col = rainbow(length(unique(all_tumours_joined_tumour_type))), # Different colors for each tumor type
        xlab = "Tumor Type", # Custom x-axis label
        ylab = "Scaled Expression Values", # Custom y-axis label
        main = "Boxplot of Scaled Expression Values by Tumor Type") # Main title


################################################################################
scale_centre <- function(x, tumour_type) {
  df <- data.frame(tumour_type = tumour_type, exp = x)
  scaled_expression <- ave(df$exp, df$tumour_type, FUN = function(y) scale(y, center = TRUE, scale = TRUE))
  return(as.numeric(scaled_expression))
}

# Apply the scaling function to each row (gene)
all_tumours_scaled_expression <- t(apply(all_tumours_joined_data_exp, 1, scale_centre, tumour_type = all_tumours_joined_tumour_type))
colnames(all_tumours_scaled_expression) <- colnames(all_tumours_joined_data_exp)



# Scaling individual tumour type data
filt.matched.ACC.exp.assay.data.scaled <- t(scale (t(filt.matched.ACC.exp.assay.data)))
filt.matched.ESCA.exp.assay.data.scaled <- t(scale (t(filt.matched.ESCA.exp.assay.data)))
filt.matched.GBM.exp.assay.data.scaled <- t(scale (t(filt.matched.GBM.exp.assay.data)))
filt.matched.LGG.exp.assay.data.scaled <- t(scale (t(filt.matched.LGG.exp.assay.data)))
filt.matched.SARC.exp.assay.data.scaled <- t(scale (t(filt.matched.SARC.exp.assay.data)))
filt.matched.BLCA.exp.assay.data.scaled <- t(scale (t(filt.matched.BLCA.exp.assay.data)))
filt.matched.KIRP.exp.assay.data.scaled <- t(scale (t(filt.matched.KIRP.exp.assay.data)))
filt.matched.PAAD.exp.assay.data.scaled <- t(scale (t(filt.matched.PAAD.exp.assay.data)))
filt.matched.TGCT.exp.assay.data.scaled <- t(scale (t(filt.matched.TGCT.exp.assay.data)))
filt.matched.LIHC.exp.assay.data.scaled <- t(scale (t(filt.matched.LIHC.exp.assay.data)))


filt.matched.ACC.exp.assay.data.scaled.df <- as.data.frame(filt.matched.ACC.exp.assay.data.scaled)
filt.matched.ESCA.exp.assay.data.scaled.df <- as.data.frame(filt.matched.ESCA.exp.assay.data.scaled)
filt.matched.GBM.exp.assay.data.scaled.df <- as.data.frame(filt.matched.GBM.exp.assay.data.scaled)
filt.matched.LGG.exp.assay.data.scaled.df <- as.data.frame(filt.matched.LGG.exp.assay.data.scaled)
filt.matched.SARC.exp.assay.data.scaled.df <- as.data.frame(filt.matched.SARC.exp.assay.data.scaled)
filt.matched.BLCA.exp.assay.data.scaled.df <- as.data.frame(filt.matched.BLCA.exp.assay.data.scaled)
filt.matched.KIRP.exp.assay.data.scaled.df <- as.data.frame(filt.matched.KIRP.exp.assay.data.scaled)
filt.matched.PAAD.exp.assay.data.scaled.df <- as.data.frame(filt.matched.PAAD.exp.assay.data.scaled)
filt.matched.TGCT.exp.assay.data.scaled.df <- as.data.frame(filt.matched.TGCT.exp.assay.data.scaled)
filt.matched.LIHC.exp.assay.data.scaled.df <- as.data.frame(filt.matched.LIHC.exp.assay.data.scaled)




filt.matched.ESCA.exp.assay.data.df.scale <- filt.matched.ESCA.exp.assay.data.df %>%
  mutate(Tumour_Type = "ESCA")
filt.matched.GBM.exp.assay.data.df.scale <- filt.matched.GBM.exp.assay.data.df %>%
  mutate(Tumour_Type = "GBM")
filt.matched.LGG.exp.assay.data.df.scale <- filt.matched.LGG.exp.assay.data.df %>%
  mutate(Tumour_Type = "LGG")


ESCA_expression_matrix <- as.matrix(filt.matched.ESCA.exp.assay.data.df.scale[, 
                                                                              -ncol(filt.matched.ESCA.exp.assay.data.df.scale)])
GBM_expression_matrix <- as.matrix(filt.matched.GBM.exp.assay.data.df.scale[, 
                                                                            -ncol(filt.matched.GBM.exp.assay.data.df.scale)])
LGG_expression_matrix <- as.matrix(filt.matched.LGG.exp.assay.data.df.scale[, 
                                                                            -ncol(filt.matched.LGG.exp.assay.data.df.scale)])

tumour_type <- filt.matched.ECSA.exp.assay.data.df.scale$Tumour_Type
tumour_type <- filt.matched.GBM.exp.assay.data.df.scale$Tumour_Type
tumour_type <- filt.matched.LGG.exp.assay.data.df.scale$Tumour_Type
