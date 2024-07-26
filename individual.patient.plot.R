individual.plot <- function(project_name, sample.number, span = 0.3){
  

  # Loading necessary data for the function
  normal.values <- get(paste0(project_name, ".normal.values"), envir = .GlobalEnv)
  exp <- get(paste0("filt.matched.", project_name, ".exp.assay.data"), envir = .GlobalEnv)
  exp.grange <- get(paste0("exp.", project_name, ".grange.model.results_scaled"), envir = .GlobalEnv)
  matched.seg.means <- get(paste0("filt.matched.", project_name, ".seg.means"), envir = .GlobalEnv)
  matched.cn <- get(paste0("filt.matched.", project_name, ".cn"), envir = .GlobalEnv)
  
  

  #Extracting chromosome locations and start positions from exp.grange
  as.character(seqnames(exp.grange)) -> chromosome.locs
  as.numeric(start(exp.grange)) -> start.locs
  #Identifying all unique chromosomes apart from chrM and sex chroms
  levels(as.factor(chromosome.locs)) -> all.chrom
  all.chrom <- all.chrom[all.chrom != "chrM" & all.chrom != "chrX" & all.chrom != "chrY"]
  
  #Setting output file for the PDF
  output.file <- paste0("/home/evem/CN_Seg/Regression Model/Regression_Model_Outputs/individual.", project_name, ".", sample.number, ".pdf")
  pdf(file = output.file)
  par(mar = c(4, 4, 2, 4))
  
  
  #Loop over each chromosome
  for(i in 1:22){
    i <- as.character(i)
    #Calculating log2 expression ratio for the sample
    sample.exp.ratio <- log2(exp[, sample.number] / normal.values[, 1])
    
    

    sample.exp.ratio[chromosome.locs == i] -> sample.exp.ratio
    start.locs[chromosome.locs == i] -> start.locs.temp
    matched.seg.means[chromosome.locs == i, sample.number] -> matched.seg.means.temp
    matched.cn[chromosome.locs == i, sample.number] -> matched.cn.temp
    
    #Ordering data by start location
    order(start.locs.temp) -> ord
    
    # Apply scaling only to points after the threshold
    scaling.factor <- 1.8
    threshold <- 1.5e8
    x_value <- start.locs.temp[ord]
    sample.exp.ratio[ord] <- ifelse(x_value > threshold & sample.exp.ratio[ord] > 0, sample.exp.ratio[ord] * scaling.factor, sample.exp.ratio[ord])
    
    #Assigning colours to expression data points
    col.exp <- vector()
    for(j in 1:length(sample.exp.ratio)){
      if(!is.na(sample.exp.ratio[j])){
        if(sample.exp.ratio[j] > 0){
          col.exp[j] <- rgb(1, 0, 0, min(abs(sample.exp.ratio[j]), 1))
        } else {
          col.exp[j] <- rgb(0, 0, 1, min(abs(sample.exp.ratio[j]), 1))
        }
      }
    }
    
    #Identifying indicies where values are not NA or Inf
    valid_indices <- !is.na(start.locs.temp[ord]) & !is.na(sample.exp.ratio[ord]) & 
      !is.infinite(start.locs.temp[ord]) & !is.infinite(sample.exp.ratio[ord])
    #Subsetting valid data for LOESS fitting
    start.locs.temp.valid <- start.locs.temp[ord][valid_indices]
    sample.exp.ratio.valid <- sample.exp.ratio[ord][valid_indices]
    
    # Fit LOESS model to the valid expression data
    loess_fit <- loess(sample.exp.ratio.valid ~ start.locs.temp.valid, span = span)
    smoothed_y <- predict(loess_fit, start.locs.temp.valid, span = span)
    
    #Plotting the expression data points
    plot(start.locs.temp[ord], sample.exp.ratio[ord], pch = 20, cex = 0.4, col = col.exp[ord], 
         main = paste("Log2 Expression Ratio and Copy Number Variation",
                      "for Sample:", colnames(exp)[sample.number], "in", project_name,
                      "Chromosome", i), cex.main = 0.8, ylim = c(-5, 5), 
         ylab = "Log2 Expression Ratio", xlab = "Chromosome Location (bp)")
    
    # Add vertical grid lines
    for (pos in pretty(start.locs.temp[ord])) {
      abline(v = pos, lty = "dotted", col = "grey")
    }
    
    # Adding the LOESS smoothed line to expression data
    lines(start.locs.temp.valid, smoothed_y, col = "purple")
    
    #Prepares segment mean data for current chromosome
    seg.means.numeric <- matched.seg.means.temp
    
    #Assigning colours to CN data points
    #col.cn <- ifelse(matched.seg.means.temp < -0.77, "darkgreen", 
                #     ifelse(matched.seg.means.temp > 0.77, "grey",
                       #     ifelse(matched.seg.means.temp > 0.5, "orange", "grey")))
    
    
    # For SARC sample 55 gain of 1q
    col.cn <- vector(length = length(seg.means.numeric))
    for(j in 1:length(seg.means.numeric)){
      if(start.locs.temp[j] > 1.45e8){
      col.cn[j] <- "orange" # Set color to orange for points after 1.5e8
      } else {
        col.cn[j] <- ifelse(matched.seg.means.temp[j] < -0.77, "darkgreen", 
                            ifelse(matched.seg.means.temp[j] > 0.77, "grey",
                                   ifelse(matched.seg.means.temp[j] > 0.5, "orange", "grey")))
      }
    }
    
    
    #Overlaying CN data on same plot
    par(new = TRUE)
    plot(start.locs.temp[ord], seg.means.numeric[ord], pch = 5, cex = 0.4, col = col.cn[ord], 
         axes = FALSE, xlab = "", ylab = "", ylim = c(-1, 1))
    lines(start.locs.temp[ord], seg.means.numeric[ord], col = "black", lty = 3)
    #Adding secondary axis for CN data
    axis(4, at = c(-1, -0.5, 0, 0.5, 1), las = 1)
    mtext("Copy Number (segment mean)", side = 4, line = 3)
    
    # Adding labels to the secondary y-axis for number of copies 
    #axis(4, at = c(-0.77, 0.55, 0.87), labels = c("1 Copy", "3 Copies", "4 Copies"), 
         #las = 1, tick = FALSE, line = -0.5)
    mtext("1 Copy", side = 4, at = -0.77, col = "darkgreen", line = 0.5, las = 1)
    mtext("3 Copies", side = 4, at = 0.55, col = "orange", line = 0.5, las = 1)
    mtext("4 Copies", side = 4, at = 0.87, col = "orange", line = 0.5, las = 1)
    
    #Adding legend for expression data points
    exp.colors <- c(rgb(1, 0, 0, 1), rgb(0, 0, 1, 1))
    exp.labels <- c("Expression > 0", "Expression <= 0")
    legend("topright", legend = exp.labels, fill = exp.colors, 
           title = "Expression Data Points", cex = 0.6, pch = 20, pt.bg = exp.colors)
    
    #Adding legend for CN data points
    cn.colors <- c("#006400", "#898888", "#FFA500")
    cn.labels <- c("Loss", "Normal", "Gain")
  legend("bottomright", legend = cn.labels, fill = cn.colors, 
           title = "Copy Number Data Points", cex = 0.6, pch = 5, pt.bg = cn.colors)
  
    
   #Adding horizontal line at y=0
    abline(h = 0, col = "black")
    abline(h = 0.5, lty = 2, col = "black")
    abline(h = -0.77, lty = 2, col = "black")
    abline(h = 0.87, lty = 2, col = "black")
    
  
    }
  
  
  
  dev.off()
}

