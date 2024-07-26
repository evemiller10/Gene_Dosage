# Define a function to create a plot for EDF
plot_edf <- function(project_name, project_grange, chromosome_of_interest) {
  
  output.file <- paste0("/home/evem/CN_Seg/Regression Model/Regression_Model_Outputs/edf_plot_", 
                        project_name, ".Chr", chromosome_of_interest, ".pdf")
  
  # Create the data frame
  plot.dataframe <- data.frame(
    chromosome.position = as.numeric(start(project_grange)), 
    genename = project_grange$gene_name,
    ensembl.id = project_grange$gene_id,
    linear.coefficient = as.numeric(project_grange$lm.seg.means.coef),
    edf = as.numeric(project_grange$gam.edf),
    pval = as.numeric(-log10(project_grange$lm.adj.p.val)),
    chromosome = as.character(seqnames(project_grange))
  )
  
  # Subset the data frame to include only the chromosome of interest
  plot.dataframe <- plot.dataframe[plot.dataframe$chromosome == chromosome_of_interest, ]
  
  # Calculate the thresholds for the top 1% highest and lowest coefficients
  threshold_high <- quantile(plot.dataframe$linear.coefficient, 0.99, na.rm = TRUE)
  threshold_low <- quantile(plot.dataframe$linear.coefficient, 0.01, na.rm = TRUE)
  
  # Adding a column to indicate whether a gene should be labelled or not
  plot.dataframe$label <- ifelse(
    (plot.dataframe$linear.coefficient >= threshold_high | plot.dataframe$linear.coefficient <= threshold_low)
    & plot.dataframe$pval > 1.3,
    plot.dataframe$genename,
    NA
  )
  
  # Filter for significant positive coefficients
  significant_positive <- plot.dataframe %>% filter(linear.coefficient > 0 & pval > 1.3)
  
  # Start pdf device
  pdf(file = output.file)
  
  # Draw scatter plot for EDF
  title_text <- paste0("Effective Degrees of Freedom (EDF) on Chromosome ", chromosome_of_interest, " in ", project_name)
  
  plot <- ggplot(plot.dataframe, aes(x = chromosome.position, y = edf)) + 
    geom_point(alpha = 0.3, color = "royalblue") + 
    geom_smooth(method = "loess", se = FALSE, color = "red", linetype = "solid",
                size = 0.5, span = 0.3) +
    geom_label_repel(aes(label = label), size = 1.5, max.overlaps = Inf) +
    scale_y_continuous(limits = c(0, max(plot.dataframe$edf, na.rm = TRUE) * 1.1)) +
    xlab("Chromosome Position (bp)") + 
    ylab("Effective Degrees of Freedom (EDF)") +
    ggtitle(title_text) + 
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.title.x = element_text(size = 10), 
      axis.title.y = element_text(size = 10),
      plot.background = element_rect(fill = "white")
    )
  
  print(plot)
  ggsave(output.file, plot)
  dev.off()
  
  return(plot)
}

# Example usage:
# Replace 'project_name', 'project_grange', and 'chromosome_of_interest' with your actual data
plot_edf("PAAD", exp.PAAD.grange.model.results_scaled, "1")  
