# Volcano plot####
EnhancedVolcano(
  limma_result,
  lab = limma_result$GeneID,
  x = 'logFC',
  y = 'P.Value',
  ylim = c(0, 6),
  title = "DSK-DS",
  subtitle = NULL,
  pCutoff = 0.05,
  pCutoffCol = 'adj.P.Val',
  FCcutoff = 0.58,
  max.overlaps = 10,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  endsConnectors = "last",
  legendPosition = "top",
  colAlpha = 0.5,
  caption = bquote( ~ Log[2] ~ "fold change cutoff: 0.58; adj.P.Val cutoff: < 0.05")
)

ggsave(
  paste0(
    DE_path,
    "Volcano_",
    "DSK-DS",
    "_",
    "Limma_norm_imputed_adj.P.Val<0.05_",
    projectname,
    "_200x200mm_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 200,
  width = 200,
  units = "mm",
  dpi = "retina"
)


# MD-Plot#####
limma_plot <- limma_result  %>%
  mutate(
    gene_type = case_when(
      `logFC` > 0.58 & `adj.P.Val` < 0.05 ~ "Up",
      `logFC` < -0.58 &
        `adj.P.Val` < 0.05 ~ "Down",
      TRUE ~ "NotSig"
    )
  )

# Create subtable of significantly regulated genes for text labeling
sig_genes <-
  filter(limma_plot, gene_type == "Up" | gene_type == "Down")


#Create vector defining coloring of data points
cols <- c("Up" = "red",
          "Down" = "#26b3ff",
          "NotSig" = "grey")

#Plot scatter
ggplot(limma_plot, aes(x = AveExpr, y = logFC)) +
  geom_point(
    aes(fill = gene_type),
    colour = "darkgrey",
    alpha = 0.8,
    shape = 21,
    size = 2
  ) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_point(
    data = sig_genes,
    # Modify shape + colour of points of interest
    colour = "black",
    alpha = 1,
    shape = 21,
    size = 2
  ) +
  geom_text_repel(
    data = sig_genes,
    # Add labels last to appear as the top layer
    aes(label = GeneID),
    size = 5,
    force = 1.5,
    nudge_y = 0,
    max.overlaps = 10,
    min.segment.length = 0
  ) +
  scale_fill_manual(values = cols) +  # Modify point colour
  labs(
    x = "Average expression",
    y = "Log2 fold change",
    fill = "Regulation"
  ) +
  theme_bw(
  )+
  theme(
    axis.text = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 16),
    plot.title = element_text(hjust = 0, face = "bold"),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 12)
  )
ggsave(
  paste0(
    DE_path,
    "MDplot_",
    "DSK-DS",
    "_",
    "Limma_norm_imputed_adj.P.Val<0.05_",
    projectname,
    "_380x200mm_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 200,
  width = 380,
  units = "mm",
  dpi = "retina"
)
