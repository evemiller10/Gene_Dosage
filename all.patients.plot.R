library(ggplot2)
library(ggrepel)


cn_effect_size_plot <- function(project_name, project_grange, chromosome_of_interest) {
 
  output.file <- paste0("/home/evem/CN_Seg/Regression Model/Regression_Model_Outputs/all.patients.", 
                        project_name, ".Chr", chromosome_of_interest, ".pdf")
  
  
   #Create the data frame
  plot.dataframe <- data.frame(
    chromosome.position = as.numeric(start(project_grange)), 
    genename = project_grange$gene_name,
    linear.coefficient = project_grange$lm.seg.means.coef,
    edf = as.numeric(project_grange$gam.edf),
    pval = -log10(project_grange$gam.adj.p.val),
    chromosome = as.character(seqnames(project_grange))
  )
  #Subset the data frame to include only the chromosome of interest
  plot.dataframe <- plot.dataframe[plot.dataframe$chromosome == chromosome_of_interest, ]
  
  # Calculate the thresholds for the top 1% highest and lowest coefficients
  threshold_high <- quantile(plot.dataframe$linear.coefficient, 0.99, na.rm = TRUE)
  threshold_low <- quantile(plot.dataframe$linear.coefficient, 0.01, na.rm = TRUE)
  
  # Adding a column to indicate whether a gene should be labelled or not
  plot.dataframe$label <- ifelse(
    plot.dataframe$linear.coefficient >= threshold_high | plot.dataframe$linear.coefficient <= threshold_low,
    plot.dataframe$genename,
    NA
  )
  
  #Cap the linear coefficients at 300 and -300
  plot.dataframe$capped.coefficient <- ifelse(plot.dataframe$linear.coefficient > 200, 200, 
                                       ifelse(plot.dataframe$linear.coefficient < -200, -200,
                                              plot.dataframe$linear.coefficient))
  #Add a column to indicate capped values
  plot.dataframe$is_capped <- ifelse(plot.dataframe$linear.coefficient > 200 |
                              plot.dataframe$linear.coefficient < -200, TRUE, FALSE)
  
  #Adding a custom transformation function
  size_of_pval <- function() {
    trans <- function(x) ifelse(x < 1.3, x, x^2)
    inv <- function(x) ifelse(x < 1.3, x, sqrt(x))
    scales::trans_new("custom", trans, inv)
  }
  
  #Start pdf device
  pdf(file = output.file)
    
  #Draw scatter plot
  title_text <- paste0("Effect of Copy Number Alterations on Gene Expression 
                       on Chromosome ", chromosome_of_interest, " in ", project_name)
  
 print(plot <- ggplot(plot.dataframe, aes(x= chromosome.position, y = capped.coefficient,
              size = pval, color = edf)) + 
              geom_point(aes(shape = is_capped), alpha = 0.4) + 
              geom_label_repel(aes(label = label), size = 2) +
              scale_x_log10() + scale_size_continuous(trans = size_of_pval(), breaks = c(1, 1.3, 2, 3),
              range = c(0.5, 3)) + scale_color_gradient(trans = "reverse") +
              xlab("Chromosome Position (bp)") + 
              ylab("Coefficient (Effect of CN on Gene Expression)") +
              labs(size = "Significance (-log10(p-value)", color = "Linearity (EDF)") 
              + geom_hline(yintercept = 0, linetype = "solid", color = "black") +
              ggtitle(title_text) + theme_minimal() + guides(shape = guide_legend(title = "Capped at y-limit")) +
              theme(legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 8),
              legend.title = element_text(size = 7)) + ylim(-200, 200))
 
 ggsave("/home/evem/CN_Seg/Regression Model/Regression_Model_Outputs/Chr1.test",
        plot)

 dev.off()
}
