#Gain of 1q
chrom1q_genes_acc <- subset(exp.ACC.grange.model.results_scaled, seqnames == "1" & start >= 1.5e8 & start <= 2.5e8)
chrom1q_genes_blca <- subset(exp.BLCA.grange.model.results_scaled, seqnames == "1" & start >= 1.5e8 & start <= 2.5e8)
chrom1q_genes_esca <- subset(exp.ESCA.grange.model.results_scaled, seqnames == "1" & start >= 1.5e8 & start <= 2.5e8)
  chrom1q_genes_gbm <- subset(exp.GBM.grange.model.results_scaled, seqnames == "1" & start >= 1.5e8 & start <= 2.5e8)
  chrom1q_genes_kirp <- subset(exp.KIRP.grange.model.results_scaled, seqnames == "1" & start >= 1.5e8 & start <= 2.5e8)
  chrom1q_genes_lihc <- subset(exp.LIHC.grange.model.results_scaled, seqnames == "1" & start >= 1.5e8 & start <= 2.5e8)
  chrom1q_genes_lgg <- subset(exp.LGG.grange.model.results_scaled, seqnames == "1" & start >= 1.5e8 & start <= 2.5e8)
  chrom1q_genes_paad <- subset(exp.PAAD.grange.model.results_scaled, seqnames == "1" & start >= 1.5e8 & start <= 2.5e8)
  chrom1q_genes_sarc <- subset(exp.SARC.grange.model.results_scaled, seqnames == "1" & start >= 1.5e8 & start <= 2.5e8)
  chrom1q_genes_tgct <- subset(exp.TGCT.grange.model.results_scaled, seqnames == "1" & start >= 1.5e8 & start <= 2.5e8)


library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Modified function to plot multiple projects
cn_effect_size_plot_scaled_multi <- function(project_names, project_granges, chromosome_of_interest) {
  # Combine data from all projects into one dataframe
  plot.dataframes <- lapply(seq_along(project_granges), function(i) {
    project_grange <- project_granges[[i]]
    project_name <- project_names[i]
    
    data.frame(
      Tumour_Type = project_name,
      chromosome.position = as.numeric(start(project_grange)), 
      genename = mcols(project_grange)$gene_name,
      ensembl.id = mcols(project_grange)$gene_id,
      linear.coefficient = as.numeric(mcols(project_grange)$lm.seg.means.coef),
      edf = as.numeric(mcols(project_grange)$gam.edf),
      pval = as.numeric(-log10(mcols(project_grange)$lm.adj.p.val)),
      chromosome = as.character(seqnames(project_grange))
    )
  }) %>% bind_rows()
  
  # Subset for the chromosome of interest
  plot.dataframes <- plot.dataframes[plot.dataframes$chromosome == chromosome_of_interest, ]

  
  # Plotting
  title_text <- paste0("Effect of Copy Number Alteration on Gene Expression on Chromosome Arm 1q ")
  
  plot <- ggplot(plot.dataframes, aes(x = chromosome.position, y = linear.coefficient, color = Tumour_Type)) + 
    geom_point(alpha = 0.3) + 
    geom_smooth(aes(group = Tumour_Type), method = "loess", se = FALSE, linetype = "solid", size = 0.8, span = 0.3) +
    xlab("Chromosome Position (bp)") + 
    ylab("Linear Coefficient (Effect of CN on Gene Expression)") +
    ggtitle(title_text) + 
    theme_minimal() + 
    scale_color_manual(values = c("olivedrab2", "magenta", "turquoise")) + # Adjust colors as needed
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_text(size = 10), 
          axis.title.y = element_text(size = 10))
  
  print(plot)
}

# Example usage
project_names <- c("GBM", "PAAD", "SARC")
project_granges <- list(chrom1q_genes_gbm, chrom1q_genes_paad, chrom1q_genes_sarc)
chromosome_of_interest <- "1"

cn_effect_size_plot_scaled_multi(project_names, project_granges, chromosome_of_interest)
