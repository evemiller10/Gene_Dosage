ggplot(fgsea_results_combined, aes(x = Tumor, y = pathway, size = size, color = color)) +
  geom_point() +
  scale_size_continuous(name = "NES", range = c(1, 5)) +
  scale_color_continuous(name = "-log10(padj)", low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Hallmark Gene Sets Pathway Analysis across Tumour Types", x = "Tumour Type", y = "Pathway")
