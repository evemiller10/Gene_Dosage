cn.effect.multiple.tumours <- function(project_names, project_granges, chromosome_of_interest) {
  plot <- ggplot()
  
  for (i in seq_along(project_names)) {
    project_name <- project_names[i]
    project_grange <- project_granges[[i]]
    
    plot.dataframe <- data.frame(
      chromosome.position = as.numeric(start(project_grange)), 
      genename = project_grange$gene_name,
      linear.coefficient = project_grange$lm.seg.means.coef,
      edf = as.numeric(project_grange$gam.edf),
      pval = -log10(project_grange$gam.adj.p.val),
      chromosome = as.character(seqnames(project_grange)),
      project = project_name)
    
    #Subset the data frame to include only the chromosome of interest
    plot.dataframe <- plot.dataframe[plot.dataframe$chromosome == chromosome_of_interest, ]
    
    plot.dataframe$capped.coefficient <- ifelse(plot.dataframe$linear.coefficient > 200, 200, 
                                                ifelse(plot.dataframe$linear.coefficient < -200, -200,
                                                       plot.dataframe$linear.coefficient))
    #Add a column to indicate capped values
    plot.dataframe$is_capped <- ifelse(plot.dataframe$linear.coefficient > 200 |
                                         plot.dataframe$linear.coefficient < -200, TRUE, FALSE)
    
    title_text <- paste0("Effect of Copy Number Alterations on Gene Expression 
                        Across Multiple Tumors for Chromosome", chromosome_of_interest)
    
    # Add the scatter plot and loess regression line to the existing plot
    plot <- plot + 
      #geom_point(data = plot.dataframe, aes(x = chromosome.position, y = capped.coefficient,
                                            #size = pval, color = edf, shape = is_capped), alpha = 0.4) +
      geom_smooth(data = plot.dataframe, aes(x = chromosome.position, y = capped.coefficient, color = project),
                  method = "loess", se = FALSE) + xlab("Chromosome Position (bp)") + 
      ylab("Coefficient (Effect of CN on Gene Expression)")  + 
      geom_hline(yintercept = 0, linetype = "solid", color = "black") + ggtitle(title_text) +
     theme_minimal() 
  }
  plot <- plot + scale_color_manual(values = 1:length(project_names), labels = project_names)
  
  
  print(plot)

  ggsave("/home/evem/CN_Seg/Regression Model/Regression_Model_Outputs/Chr1.test.pdf",
         plot)
}
