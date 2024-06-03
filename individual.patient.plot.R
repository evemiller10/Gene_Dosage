library(ggplot2)
library(scales)


individual.plot<-function(project_name, sample.number, span = 0.05){
  
  #Loading in necessary data for the function
  normal.values <- get(paste0(project_name, ".normal.values"), envir = .GlobalEnv)
  exp <- get(paste0("filt.matched.", project_name, ".exp.assay.data"), envir = .GlobalEnv)
  exp.grange <- get(paste0("filt.exp.", project_name, ".grange.2"), envir = .GlobalEnv)
  matched.seg.means <- get(paste0("filt.matched.", project_name, ".seg.means"), envir = .GlobalEnv)
  matched.cn <- get(paste0("filt.matched.", project_name, ".cn"), envir = .GlobalEnv)
  
  
  as.character(seqnames(exp.grange))->chromosome.locs
  as.numeric(start(exp.grange))->start.locs
  levels(as.factor(chromosome.locs))->all.chrom
  all.chrom <- all.chrom[all.chrom!="chrM"&all.chrom!="chrX"&all.chrom!="chrY"]
  output.file <- paste0("/home/evem/CN_Seg/Regression Model/Regression_Model_Outputs/individual.", project_name, ".", sample.number, ".pdf")
  pdf(file = output.file)
  par(mar = c(4, 4, 2, 4))
  
  for(i in c("1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")){
    log(exp[,sample.number]/normal.values[,1],2) -> sample.exp.ratio
    sample.exp.ratio[chromosome.locs==i] -> sample.exp.ratio
    start.locs[chromosome.locs==i] -> start.locs.temp
    matched.seg.means[chromosome.locs==i,sample.number]->matched.seg.means.temp
    matched.cn[chromosome.locs==i,sample.number]->matched.cn.temp
    
    col.exp <- vector()
    for(j in 1:length(sample.exp.ratio)){
      if(!is.na(sample.exp.ratio[j])){
        if(sample.exp.ratio[j] > 0){
          # For positive values, use red color (rgb(1,0,0,alpha))
          # Alpha is determined by the value of sample.exp.ratio[j], with a maximum of 1
          col.exp[j] <- rgb(1, 0, 0, min(abs(sample.exp.ratio[j]), 1))
        } else {
          # For negative values, use blue color (rgb(0,0,1,alpha))
          # Alpha is determined by the value of sample.exp.ratio[j], with a maximum of 1
          col.exp[j] <- rgb(0, 0, 1, min(abs(sample.exp.ratio[j]), 1))
        }
      }
    }

    plot(start.locs.temp, sample.exp.ratio, pch = 20, cex = 0.4, col = col.exp, 
         main = paste("Log2 Expression Ratio and Copy Number Variation",
                      "for Sample:", colnames(exp)[sample.number], "in", project_name,
                      "Chromosome", i), cex.main = 0.8, ylim = c(-5, 5), 
         ylab = "Log2 Expression Ratio", xlab = "Chromosome Location (bp)")
   
    # Add grid lines
    grid(lty = "dotted", col = "black", lwd = 0.5)
    
    col.cn<-vector()
    for(j in 1:length(sample.exp.ratio)){
      col.cn[j] <- ifelse(matched.seg.means.temp[j] < -0.3,"#2739D7", 
                          ifelse(matched.seg.means.temp[j] <= 0.583,"#898888", 
                          ifelse(matched.seg.means.temp[j] <= 1,"#FFA500", "#FFC0CB")))
    }
    
    seg.means.numeric <- matched.seg.means.temp * 10
    
    par(new = TRUE)
    plot(start.locs.temp, seg.means.numeric, pch = 5, cex = 0.4, col = col.cn, 
         axes = FALSE, xlab = "", ylab = "", ylim = c(-10, 10))
    axis(4, at = c(-10, -5, 0, 5, 10), las = 1)
    mtext("Copy Number (segment mean)", side = 4, line = 3)
    
    points(start.locs.temp, matched.seg.means.temp*10, pch = 20, cex = 0.2, col = col.cn)
    
    # Define the colors and labels for the expression data points
    exp.colors <- c(rgb(1, 0, 0, 1), rgb(0, 0, 1, 1))
    exp.labels <- c("Expression > 0", "Expression <= 0")
    
    # Add the legend for the expression data points
    legend("topright", legend = exp.labels, fill = exp.colors, 
           title = "Expression Data Points", cex = 0.6, pch = 20, pt.bg = exp.colors)
    
    # Define the colors and labels for the copy number data points
    cn.colors <- c("#2739D7", "#898888", "#FFA500", "#FFC0CB")
    cn.labels <- c("Loss", "Normal", "Gain", "Amplification")
    
    # Add the legend for the copy number data points
    legend("bottomright", legend = cn.labels, fill = cn.colors, 
           title = "Copy Number Data Points", cex = 0.6, pch = 5, pt.bg = cn.colors)
    
    
    abline(h = 0)
    
  }
  dev.off()
}

###############################################################################
#Trying to make plot with ggplot but graph using basic R plotting is easier to
#interpret 

individual.ggplot<-function(project_name, sample.number, span = 0.05){
  
  #Loading in necessary data for the function
  normal.values <- get(paste0(project_name, ".normal.values"), envir = .GlobalEnv)
  exp <- get(paste0("filt.matched.", project_name, ".exp.assay.data"), envir = .GlobalEnv)
  exp.grange <- get(paste0("filt.exp.", project_name, ".grange.2"), envir = .GlobalEnv)
  matched.seg.means <- get(paste0("filt.matched.", project_name, ".seg.means"), envir = .GlobalEnv)
  matched.cn <- get(paste0("filt.matched.", project_name, ".cn"), envir = .GlobalEnv)
  
  
  as.character(seqnames(exp.grange))->chromosome.locs
  as.numeric(start(exp.grange))->start.locs
  levels(as.factor(chromosome.locs))->all.chrom
  all.chrom <- all.chrom[all.chrom!="chrM"&all.chrom!="chrX"&all.chrom!="chrY"]
  output.file <- paste0("/home/evem/CN_Seg/Regression Model/Regression_Model_Outputs/individual.", project_name, ".", sample.number, ".pdf")
  pdf(file = output.file)
  par(mar = c(4, 4, 2, 4))
 
  
  for(i in c("1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")){
    log(exp[,sample.number]/normal.values[,1],2) -> sample.exp.ratio
    sample.exp.ratio[chromosome.locs==i] -> sample.exp.ratio
    start.locs[chromosome.locs==i] -> start.locs.temp
    matched.seg.means[chromosome.locs==i,sample.number]->matched.seg.means.temp
    matched.cn[chromosome.locs==i,sample.number]->matched.cn.temp
    seg.means.numeric <- matched.seg.means.temp * 10

    col.exp<-vector()
    for(j in 1:length(sample.exp.ratio)){
      col.exp[j] <- ifelse(sample.exp.ratio[j]>1|sample.exp.ratio[j]< -1,rgb(1,0,0,1), 
                           ifelse(sample.exp.ratio[j]>0.5|sample.exp.ratio[j]< -0.5,rgb(1,0,0,0.3), 
                                  rgb(1,0,0,0.1)))
    }

    col.cn<-vector()
    for(j in 1:length(sample.exp.ratio)){
      col.cn[j] <- ifelse(matched.seg.means.temp[j] < -0.3,"#2739D7", 
                          ifelse(matched.seg.means.temp[j] <= 0.583,"#898888", 
                                 ifelse(matched.seg.means.temp[j] <= 1,"#FFA500", "#FFC0CB")))
    }
    
    # Create a data frame for ggplot
    df <- data.frame(
      start.locs.temp = start.locs.temp,
      sample.exp.ratio = sample.exp.ratio,
      seg.means.numeric = seg.means.numeric,
      col.exp = col.exp,
      col.cn = col.cn
    )
    
    # Create the plot
    p <- ggplot(df, aes(x = start.locs.temp)) +
      geom_point(aes(y = sample.exp.ratio, color = col.exp), size = 0.2) +
      geom_point(aes(y = seg.means.numeric, color = col.cn), size = 0.2) +
      scale_color_identity() +
      labs(
        title = paste("Log2 Expression Ratio and Copy Number Variation",
                      "for Sample:", colnames(exp)[sample.number], "in", project_name,
                      "Chromosome", i),
        x = "Chromosome Location (bp)",
        y = "Log2 Expression Ratio"
      ) +
      theme(plot.title = element_text(size = 8)) +
      ylim(-10, 10) +
      scale_y_continuous(sec.axis = sec_axis(~ . * 10, name = "Copy Number (segment mean)"))
    
    # Save the plot
    ggsave(paste0("/home/evem/CN_Seg/Regression Model/Regression_Model_Outputs/individual.", 
                  project_name, ".", sample.number, ".chromosome", i, ".pdf"), p)
  }
    return(p)
}
