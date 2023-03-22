### Get data from exported excel####

### Filter only for relevant abundances, omit DS ####
mydata_noNaN <- psite_table %>%
  select(contains("imputed_intensity")) %>%
  select(contains(c("DSK", "SK", "K")))

### PCA####
#1. Prepare data for PCA
#Remove NaN
mydata_noNaN <- na.omit(dat_norm_impute)

#Transpose data
mydata_noNaN <- t(mydata_noNaN)

#2. Perform PCA
mydata_pca <- prcomp(mydata_noNaN)
summ <- summary(mydata_pca)

#3. Create new dataframe containing only PC1 and PC2
mydata_PC <- as.data.frame(mydata_pca$x)

#4. Add new column with group labels
mydata_PC$Group <- rownames(mydata_PC)

#5. Split group label column to Sample_type and Replicate
mydata_PC <-
  separate(
    data = mydata_PC,
    col = "Group",
    into = c(1,2,3,4,"Sample_type", "Label", "Replicate"),
    sep = "_"
  )

#6. Change Sample_type to factor
mydata_PC <-
  mutate(.data = mydata_PC,
         Sample_type = factor(Sample_type, levels = c("DSK",
                                                      "DS",
                                                      "SK",
                                                      "K")))


#7. Calculate percent variance values for x and y labs, round Percentage to 1 digit
percentage_var <-
  round(mydata_pca$sdev ^ 2 / sum(mydata_pca$sdev ^ 2) * 100, 1)
percentage_var <-
  paste(colnames(mydata_PC),
        "(",
        paste(as.character(percentage_var), "%", ")", sep = ""))


#8. Plot data with ggplot2
ggplot(data = mydata_PC, aes(x = PC1, y = PC2, color = Sample_type)) +
  geom_point(size = 5, show.legend = TRUE) +
  geom_text_repel(
    aes(label = paste0(Sample_type, "_", Replicate)),
    show.legend = FALSE,
    size = 5,
    fontface = "bold",
    force = 0.5,
    max.overlaps = Inf,
    min.segment.length = 0,
    point.padding = 1,
    box.padding = 1,
    segment.linetype = "dashed"
  ) +
  scale_color_manual(values = heatcolors2[[1]],
                     name = "Construct") +
  theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank()
  ) +
  xlab(percentage_var[1]) + ylab(percentage_var[2])

#9. Save plot to file
ggsave(
  paste0(QC_path,
         "P-sites_PCA_",
         "enriched_norm_log2_intensities+impute_L_woDS",
         "_200x165mm_",
         projectname,
         "_",
         format(Sys.time(), '%Y%m%d_%H%M%S'),
         ".png"
  ),
  height = 165,
  width = 200,
  units = "mm",
  dpi = "retina"
)

ggsave(
  paste0(QC_path,
         "P-sites_PCA_",
         "enriched_norm_log2_intensities+impute_L_woDS",
         "_200x165mm_",
         projectname,
         "_",
         format(Sys.time(), '%Y%m%d_%H%M%S'),
         ".svg"
  ),
  height = 165,
  width = 200,
  units = "mm",
  dpi = "retina")