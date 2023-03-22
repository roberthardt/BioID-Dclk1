#Cleaning of workspace, console and plot area

# Clear plots
if (!is.null(dev.list()))
  dev.off()
# Clear console
cat("\014")
# Clean workspace
rm(list = ls())
#-------------------------------------------------------------------------------
#Define needed package
package.names <-
  c(
    "tidyverse",
    "DEqMS",
    "tcltk",
    "matrixStats",
    "scales",
    "svglite",
    "ggsci",
    "ggrepel",
    "BiocManager",
    "pheatmap",
    "heatmaply",
    "openxlsx",
    "BioVenn",
    "reshape2",
    "RColorBrewer"
  )

#Install/Load packages

is_installed <-
  function(mypkg)
    is.element(mypkg, installed.packages()[, 1])

load_or_install <- function(package_names)
{
  for (package_name in package_names)
  {
    if (!is_installed(package_name))
    {
      install.packages(package_name, repos = "http://lib.stat.cmu.edu/R/CRAN")
    }
    library(
      package_name,
      character.only = TRUE,
      quietly = TRUE,
      verbose = FALSE
    )
  }
}

load_or_install(package.names)

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

# Function for opening a pop-up window (INDEPENDENT of operating system)
# for selecting the working directory
choose_directory = function(caption = 'Select data directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = "Select folder")
  } else {
    tk_choose.dir(caption = "Select folder")
  }
}


# Function for opening a pop-up window (INDEPENDENT of operating system)
# for selecting a file
choose_file = function(caption = 'Select file') {
  if (exists('utils::choose.files')) {
    choose.files(caption = "Select files")
  } else {
    tk_choose.files(caption = "Select files")
  }
}

#Set path of working directory
setwd(choose_directory())

#Select file
protein <-
  read.table(
    file = choose_file(),
    stringsAsFactors = FALSE,
    header = TRUE,
    quote = "",
    comment.char = "",
    sep = "\t"
  )


#Set project name. This name is added to all result tables and plots. 
#Please change for every project/dataset to be analyzed.
projectname <- "03_BioID-DCLK_proteome"

# remove reverse, contaminant and identified by site entries
p_filter <- protein %>%
  filter(Reverse != "+") %>%
  filter(Potential.contaminant != "+") %>%
  filter(Only.identified.by.site != "+")


# Assign majority protein id as row names
rownames(p_filter) <- p_filter$Majority.protein.IDs

# Filter for expression columns
dat <- select(p_filter, contains(c(
  "Intensity.L.BioID",
  "Intensity.H.BioID"
)))

# Set 0 values to NA
dat[dat == 0] <- NA

#rename columns
colnames(dat) <- c(
  "DS_L_01",
  "DS_L_02",
  "DS_L_03",
  "DSK_L_01",
  "DSK_L_02",
  "DSK_L_03",
  "K_L_01",
  "K_L_02",
  "K_L_03",
  "SK_L_01",
  "SK_L_02",
  "SK_L_03",
  "DS_H_01",
  "DS_H_02",
  "DS_H_03",
  "DSK_H_01",
  "DSK_H_02",
  "DSK_H_03",
  "K_H_01",
  "K_H_02",
  "K_H_03",
  "SK_H_01",
  "SK_H_02",
  "SK_H_03"
)

#Log2 transform data
dat_log <- dat %>%
  dplyr::mutate(across(.cols = everything(),
                       .fns = log2))


# Generate some basic qc plots (box, density) in ggplot. 
# Therefore convert to long format first
dat_log_long <- dat_log %>%
  pivot_longer(
    cols = everything(),
    names_to = c("sample", "label", "replicate"),
    names_sep = "\\_",
    values_to = "log2_intensity"
  ) %>%
  unite(col = "sampleID",
        c(sample, replicate),
        remove = FALSE) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  mutate(label = factor(label, levels = unique(label))) %>%
  mutate(replicate = factor(replicate, levels = unique(replicate))) %>%
  mutate(sampleID = factor(sampleID, levels = unique(sampleID)))


#reorder factor levels to achieve different sorting in plots
levels(dat_log_long$sample) <- c("DSK", "DS", "SK", "K")
levels(dat_log_long$replicate) <- c("01", "02", "03")
levels(dat_log_long$label) <- c("L", "H")
levels(dat_log_long$sampleID) <- c(
  "DSK_01",
  "DSK_02",
  "DSK_03",
  "DS_01",
  "DS_02",
  "DS_03",
  "SK_01",
  "SK_02",
  "SK_03",
  "K_01",
  "K_02",
  "K_03"
)

#Make plots

#Boxplot

# Calc mean log2 intensities
stat_log2_inten <- dat_log_long %>%
  group_by(label) %>%
  summarize(
    mean = mean(log2_intensity, na.rm = TRUE),
    median = median(log2_intensity, na.rm = TRUE)
  )


# Plot with ggplot
ggplot(dat_log_long, aes(x = sampleID,
                         y = log2_intensity,
                         fill = sample)) +
  geom_boxplot() +
  facet_wrap(vars(label)) +
  geom_hline(data = stat_log2_inten,
             aes(yintercept = median,
                 color = "median"),
             linetype = "dashed") +
  labs(x = "Sample type",
       y = "log2 intensity") +
  scale_fill_brewer(palette = "Dark2",
                    name = "Construct") +
  scale_color_manual(name = "Statistiscs",
                     values = "red") +
  theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(
      face = "bold",
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

#Save boxplot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(
    "Boxplot_log2_intensities_L+H_400x150_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 150,
  width = 400,
  units = "mm",
  dpi = "retina"
)

ggsave(
  paste0(
    "Boxplot_log2_intensities_L+H_400x150_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".svg"
  ),
  height = 150,
  width = 400,
  units = "mm",
  dpi = "retina"
)


#Density plots
# Calc mean log2 intensities
stat_log2_inten <- dat_log_long %>%
  group_by(label, sampleID) %>%
  summarize(
    mean = mean(log2_intensity, na.rm = TRUE),
    median = median(log2_intensity, na.rm = TRUE)
  )

#Plot with ggplot
ggplot(dat_log_long, aes(x = log2_intensity,
                         fill = sample)) +
  geom_density(alpha = 0.75) +
  facet_grid(label ~ sampleID) +
  geom_vline(data = stat_log2_inten,
             aes(xintercept = median,
                 color = "median"),
             linetype = "dashed") +
  labs(x = "log2 intensity",
       y = "density") +
  scale_fill_brewer(palette = "Dark2",
                    name = "Construct") +
  scale_color_manual(name = "Statistiscs",
                     values = "red") +
  theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(
      face = "bold",
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

#Save boxplot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(
    "Density_log2_intensities_L+H_400x100_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 100,
  width = 400,
  units = "mm",
  dpi = "retina"
)

ggsave(
  paste0(
    "Density_log2_intensities_L+H_400x100_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".svg"
  ),
  height = 100,
  width = 400,
  units = "mm",
  dpi = "retina"
)


#Create MA plots against global average
A <- rowMeans(dat_log, na.rm = TRUE)
M <- sweep(dat_log, 1, rowMeans(dat_log, na.rm = TRUE))

pdf(
  file = paste0(
    "MA-plots_nonorm",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".pdf"
  ),
  paper = "a4r",
  width = 11
)
par(mfrow = c(2, 3))
for (j in colnames(dat_log)) {
  plot(
    x = A,
    y = M[, j],
    main = j,
    xlab = "Average log-expression",
    ylab = "Expression log-ratio (this sample vs av. others)",
    pch = 20
  )
  lines(predict(loess(M[, j] ~ A)), col = "red")
  abline(h = 0, col = "black")
}
dev.off()

#Count number of NaN per condition.
dat_log$na_count_DS_L = apply(dat_log, 1, function(x)
  sum(is.na(x[c(1:3)])))
dat_log$na_count_DSK_L = apply(dat_log, 1, function(x)
  sum(is.na(x[c(4:6)])))
dat_log$na_count_K_L = apply(dat_log, 1, function(x)
  sum(is.na(x[c(7:9)])))
dat_log$na_count_SK_L = apply(dat_log, 1, function(x)
  sum(is.na(x[c(10:12)])))


# Filter protein table for valid values. Here we want a minimum of 2 values 
# in any Light sample group.
dat_log_filter <- dat_log %>%
  filter(if_any(contains("na_count"),
                ~ . <= 1)) %>%
  select(c(1:24))

# Impute missing heavy signals by normal distribution
#First store position info of imputed values in a separate data frame
na_pos <- is.na(dat_log_filter) %>%
  as.data.frame()

#Perform imputation
dat_log_impute <- dat_log_filter %>%
  as.matrix() %>%
  MsCoreUtils::impute_matrix(method = "MinProb") %>%
  as.data.frame()

# Plot density plots for data after imputation
# Generate some basic qc plots (box, density) in ggplot. 
# Therefore convert to long format first
dat_log_impute_long <- dat_log_impute %>%
  pivot_longer(
    cols = everything(),
    names_to = c("sample", "label", "replicate"),
    names_sep = "\\_",
    values_to = "log2_intensity"
  ) %>%
  unite(col = "sampleID",
        c(sample, replicate),
        remove = FALSE) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  mutate(label = factor(label, levels = unique(label))) %>%
  mutate(replicate = factor(replicate, levels = unique(replicate))) %>%
  mutate(sampleID = factor(sampleID, levels = unique(sampleID)))

#Add ID column
dat_log_impute_long <- mutate(.data = dat_log_impute_long,
                              id = row.names(dat_log_impute_long))

#Pivot also the NA-Positions to the long format
na_pos_long <- na_pos %>%
  pivot_longer(
    cols = everything(),
    names_to = c("sample", "label", "replicate"),
    names_sep = "\\_",
    values_to = "imputed"
  ) %>%
  unite(col = "sampleID",
        c(sample, replicate),
        remove = FALSE) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  mutate(label = factor(label, levels = unique(label))) %>%
  mutate(replicate = factor(replicate, levels = unique(replicate))) %>%
  mutate(sampleID = factor(sampleID, levels = unique(sampleID))) %>%
  mutate(imputed = factor(imputed, levels = unique(imputed)))

#Add ID column
na_pos_long <- mutate(.data = na_pos_long,
                      id = row.names(na_pos_long))

#Merge expression data frame and na_position frame
dat_log_impute_merge <- merge(x = dat_log_impute_long,
                              y = na_pos_long[, c("id", "imputed")],
                              by = "id")

#reorder factor levels to achieve different sorting in plots
levels(dat_log_impute_merge$sample) <- c("DSK", "DS", "SK", "K")
levels(dat_log_impute_merge$replicate) <- c("01", "02", "03")
levels(dat_log_impute_merge$label) <- c("L", "H")

levels(dat_log_impute_merge$sampleID) <-
  c(
    "DSK_01",
    "DSK_02",
    "DSK_03",
    "DS_01",
    "DS_02",
    "DS_03",
    "SK_01",
    "SK_02",
    "SK_03",
    "K_01",
    "K_02",
    "K_03"
  )


#Make plots

#Boxplot

# Calc mean log2 intensities
stat_log2_inten <- dat_log_impute_merge %>%
  group_by(label) %>%
  summarize(
    mean = mean(log2_intensity, na.rm = TRUE),
    median = median(log2_intensity, na.rm = TRUE)
  )

# Plot with ggplot
ggplot(dat_log_impute_merge,
       aes(x = sampleID,
           y = log2_intensity,
           fill = sample)) +
  geom_boxplot() +
  facet_wrap(vars(label)) +
  geom_hline(data = stat_log2_inten,
             aes(yintercept = median,
                 color = "median"),
             linetype = "dashed") +
  labs(x = "Sample type",
       y = "log2 intensity") +
  scale_fill_brewer(palette = "Dark2",
                    name = "Construct") +
  scale_color_manual(name = "Statistiscs",
                     values = "red") +
  theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(
      face = "bold",
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

#Save boxplot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(
    "Boxplot_log2_intensities+impute_L+H_400x150_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 150,
  width = 400,
  units = "mm",
  dpi = "retina"
)

ggsave(
  paste0(
    "Boxplot_log2_intensities+impute_L+H_400x150_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".svg"
  ),
  height = 150,
  width = 400,
  units = "mm",
  dpi = "retina"
)


#Density plots
# Calc mean log2 intensities
stat_log2_inten <- dat_log_impute_merge %>%
  group_by(label, sampleID) %>%
  summarize(
    mean = mean(log2_intensity, na.rm = TRUE),
    median = median(log2_intensity, na.rm = TRUE)
  )

#Plot with ggplot
ggplot(dat_log_impute_merge, aes(x = log2_intensity,
                                 fill = imputed)) +
  geom_density(alpha = 0.75) +
  facet_grid(label ~ sampleID) +
  geom_vline(data = stat_log2_inten,
             aes(xintercept = median,
                 color = "median"),
             linetype = "dashed") +
  labs(x = "log2 intensity",
       y = "density") +
  scale_fill_brewer(palette = "Dark2",
                    name = "Imputed",
                    direction = -1) +
  scale_color_manual(name = "Statistiscs",
                     values = "red") +
  theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(
      face = "bold",
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

#Save boxplot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(
    "Density_log2_intensities+impute_L+H_400x100_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 100,
  width = 400,
  units = "mm",
  dpi = "retina"
)

ggsave(
  paste0(
    "Density_log2_intensities+impute_L+H_400x100_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".svg"
  ),
  height = 100,
  width = 400,
  units = "mm",
  dpi = "retina"
)

## Generate L/H ratios
df_ratios <- dat_log_impute[, 1:12] - dat_log_impute[, 13:24]

## Count number of L/H ratios >1 per condition.
df_ratios$upreg_count_DS <- apply(df_ratios,
                                  1,
                                  function(x)
                                    sum(x[c(1:3)] > 1))
df_ratios$upreg_count_DSK <- apply(df_ratios,
                                   1,
                                   function(x)
                                     sum(x[c(4:6)] > 1))
df_ratios$upreg_count_K <- apply(df_ratios,
                                 1,
                                 function(x)
                                   sum(x[c(7:9)] > 1))
df_ratios$upreg_count_SK <- apply(df_ratios,
                                  1,
                                  function(x)
                                    sum(x[c(10:12)] > 1))

## Filter for enriched proteins. (at least two log2 ratios > 1)
df_ratios_filter <- df_ratios %>%
  filter(if_any(contains("upreg_count"),
                ~ . >= 2))
## Generate vector with protein accessions of enriched proteins
enriched_ID <- row.names(df_ratios_filter)

#Use rownames from filtered df_ratios to fetch the corresponding light log2 intensities
dat_log_enriched <- dat_log_filter[enriched_ID,
                                   c(1:12)]


##Plot light log2 intensities in boxplot and density plots of the enriched proteins

## Convert to long format
dat_log_enriched_long <- dat_log_enriched %>%
  pivot_longer(
    cols = everything(),
    names_to = c("sample", "label", "replicate"),
    names_sep = "\\_",
    values_to = "log2_intensity"
  ) %>%
  unite(col = "sampleID",
        c(sample, replicate),
        remove = FALSE) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  mutate(replicate = factor(replicate, levels = unique(replicate))) %>%
  mutate(sampleID = factor(sampleID, levels = unique(sampleID)))

#reorder factor levels to achieve different sorting in plots
levels(dat_log_enriched_long$sample) <- c("DSK", "DS", "SK", "K")
levels(dat_log_enriched_long$replicate) <- c("01", "02", "03")
levels(dat_log_enriched_long$sampleID) <-
  c(
    "DSK_01",
    "DSK_02",
    "DSK_03",
    "DS_01",
    "DS_02",
    "DS_03",
    "SK_01",
    "SK_02",
    "SK_03",
    "K_01",
    "K_02",
    "K_03"
  )
#Boxplot

# Calc mean log2 intensities
stat_log2_inten <- dat_log_enriched_long %>%
  summarize(
    mean = mean(log2_intensity, na.rm = TRUE),
    median = median(log2_intensity, na.rm = TRUE)
  )


# Plot with ggplot
ggplot(dat_log_enriched_long,
       aes(x = sampleID,
           y = log2_intensity,
           fill = sample)) +
  geom_boxplot() +
  geom_hline(data = stat_log2_inten,
             aes(yintercept = median,
                 color = "median"),
             linetype = "dashed") +
  labs(x = "Sample type",
       y = "log2 intensity") +
  scale_fill_brewer(palette = "Dark2",
                    name = "Construct") +
  scale_color_manual(name = "Statistiscs",
                     values = "red") +
  theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(
      face = "bold",
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

#Save boxplot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(
    "Boxplot_enriched_log2_intensities_L_400x150_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 150,
  width = 400,
  units = "mm",
  dpi = "retina"
)

ggsave(
  paste0(
    "Boxplot_enriched_log2_intensities_L_400x150_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".svg"
  ),
  height = 150,
  width = 400,
  units = "mm",
  dpi = "retina"
)


#Density plots
#Plot with ggplot
ggplot(dat_log_enriched_long, aes(x = log2_intensity,
                                  color = sampleID)) +
  geom_density(alpha = 0.75) +
  labs(x = "log2 intensity",
       y = "density") +
  scale_color_manual(values = rep(
    brewer.pal(n = 4,
               name = "Dark2"),
    times = 2,
    each = 3
  )) +
  theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold",
                              size = 14),
    axis.text.x = element_text(
      face = "bold",
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

#Save density plot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(
    "Density_enriched_log2_intensities_L_180x120_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 120,
  width = 180,
  units = "mm",
  dpi = "retina"
)

ggsave(
  paste0(
    "Density_enriched_log2_intensities_L_180x120_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".svg"
  ),
  height = 120,
  width = 180,
  units = "mm",
  dpi = "retina"
)

##----------------------------------------------
## Check out which normalization method to use best with NormalizerDE
library(NormalyzerDE)

#Generate design table & write to file
design_MA <- data.frame(
  sample = c(
    "DS_L_01",
    "DS_L_02",
    "DS_L_03",
    "DSK_L_01",
    "DSK_L_02",
    "DSK_L_03",
    "K_L_01",
    "K_L_02",
    "K_L_03",
    "SK_L_01",
    "SK_L_02",
    "SK_L_03"
  ),
  group = c(rep(c(
    "DS", "DSK", "K", "SK"
  ),
  each = 3)),
  batch = c(rep(c(1:3),
                times = 4))
) %>%
  write_tsv(file = "designNormalyzerDE.tsv")

##Generate input data & write to file.
# Input data has to be de-logarithmized before
input_data <- 2 ^ dat_log_enriched  %>%
  mutate(ProtID = row.names(dat_log_enriched),
         .before = 1) %>%
  write_tsv(file = "inputNormalyzerDE.tsv")

##Run normalyzerDE
NormalyzerDE::normalyzer(
  jobName = "normalizerDE_eval",
  designPath = "designNormalyzerDE.tsv",
  dataPath = "inputNormalyzerDE.tsv",
  outputDir = getwd()
)

## Result: ==> basically every method is better then simple log2. But RLR,
## cyclicLoess and vsn perform best. I opt for CyclicLoess.

##_________________________________________________________________
## Perform normalization
dat_norm <-
  normalizeBetweenArrays(dat_log_enriched, method = "cyclicloess")
##-----------------------------------------------------


## Plot normalized data in Boxplot, Density plot and PCA
## Convert to long format
dat_norm_long <- dat_norm %>%
  as.data.frame() %>%
  mutate(ProtID = row.names(dat_norm)) %>%
  pivot_longer(
    cols = c(1:12),
    names_to = c("sample", "label", "replicate"),
    names_sep = "\\_",
    values_to = "log2_intensity"
  ) %>%
  unite(col = "sampleID",
        c(sample, replicate),
        remove = FALSE) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  mutate(replicate = factor(replicate, levels = unique(replicate))) %>%
  mutate(sampleID = factor(sampleID, levels = unique(sampleID)))

#reorder factor levels to achieve different sorting in plots
levels(dat_norm_long$sample) <- c("DSK", "DS", "SK", "K")
levels(dat_norm_long$replicate) <- c("01", "02", "03")
levels(dat_norm_long$sampleID) <- c(
  "DSK_01",
  "DSK_02",
  "DSK_03",
  "DS_01",
  "DS_02",
  "DS_03",
  "SK_01",
  "SK_02",
  "SK_03",
  "K_01",
  "K_02",
  "K_03"
)
#Boxplot

# Calc mean log2 intensities
stat_log2_inten <- dat_norm_long %>%
  summarize(
    mean = mean(log2_intensity, na.rm = TRUE),
    median = median(log2_intensity, na.rm = TRUE)
  )


# Plot with ggplot
ggplot(dat_norm_long, aes(x = sampleID,
                          y = log2_intensity,
                          fill = sample)) +
  geom_boxplot() +
  geom_hline(data = stat_log2_inten,
             aes(yintercept = median,
                 color = "median"),
             linetype = "dashed") +
  labs(x = "Sample type",
       y = "log2 intensity") +
  scale_fill_brewer(palette = "Dark2",
                    name = "Construct") +
  scale_color_manual(name = "Statistiscs",
                     values = "red") +
  theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(
      face = "bold",
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

#Save boxplot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(
    "Boxplot_enriched_norm_log2_intensities_L_400x150_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 150,
  width = 400,
  units = "mm",
  dpi = "retina"
)

ggsave(
  paste0(
    "Boxplot_enriched_norm_log2_intensities_L_400x150_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".svg"
  ),
  height = 150,
  width = 400,
  units = "mm",
  dpi = "retina"
)


#Density plots
#Plot with ggplot
ggplot(dat_norm_long, aes(x = log2_intensity,
                          color = sampleID)) +
  geom_density(alpha = 0.75) +
  labs(x = "log2 intensity",
       y = "density") +
  scale_color_manual(values = rep(
    brewer.pal(n = 4,
               name = "Dark2"),
    times = 2,
    each = 3
  )) +
  theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold",
                              size = 14),
    axis.text.x = element_text(
      face = "bold",
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

#Save density plot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(
    "Density_enriched_norm_log2_intensities_L_180x120_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 120,
  width = 180,
  units = "mm",
  dpi = "retina"
)

ggsave(
  paste0(
    "Density_enriched_norm_log2_intensities_L_180x120_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".svg"
  ),
  height = 120,
  width = 180,
  units = "mm",
  dpi = "retina"
)

########### PCA


#1. Prepare data for PCA
#Remove NaN
mydata_noNaN <- na.omit(dat_norm)

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
    into = c("Sample_type", "Label", "Replicate"),
    sep = "_"
  )

#6. Change Sample_type to factor
mydata_PC <-
  mutate(.data = mydata_PC,
         Sample_type = factor(Sample_type, levels = unique(Sample_type)))


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
  scale_color_brewer(palette = "Dark2",
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
  paste0(
    "PCA_",
    "enriched_norm_log2_intensities_L",
    projectname,
    "_200x165mm_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 165,
  width = 200,
  units = "mm",
  dpi = "retina"
)

ggsave(
  paste0(
    "PCA_",
    "enriched_norm_log2_intensities_L",
    projectname,
    "_200x165mm_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".svg"
  ),
  height = 165,
  width = 200,
  units = "mm",
  dpi = "retina"
)

##Plot myc-BioID in barplot
BioID <- dat_norm_long %>% filter(ProtID == "myc-BioID")

ggplot(BioID, aes(x = sampleID,
                  y = log2_intensity,
                  fill = sample)) +
  geom_col() +
  scale_fill_brewer(palette = "Dark2",
                    name = "Construct") +
  geom_text(aes(label = signif(log2_intensity,
                               digits = 3)),
            vjust = 2)

ggsave(
  paste0(
    "Barplot_",
    "myc-BioID_norm_log2_intensities_L",
    projectname,
    "_75x100mm_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 75,
  width = 150,
  units = "mm",
  dpi = "retina"
)
##---------------------------------------------------------------
##Find proteins which are exclusively quantified in one sample group
#and store in a separate table

DS_exclusive <- dat_norm %>%
  as.data.frame() %>%
  filter(if_all(.cols = c(4:12),
                .fns = ~ is.na(.x))) %>%
  select(starts_with("DS_L"))

DSK_exclusive <- dat_norm %>%
  as.data.frame() %>%
  filter(if_all(.cols = c(1:3, 7:12),
                .fns = ~ is.na(.x))) %>%
  select(starts_with("DSK_L"))

K_exclusive <- dat_norm %>%
  as.data.frame() %>%
  filter(if_all(.cols = c(1:6, 10:12),
                .fns = ~ is.na(.x))) %>%
  select(starts_with("K_L"))

SK_exclusive <- dat_norm %>%
  as.data.frame() %>%
  filter(if_all(.cols = c(1:9),
                .fns = ~ is.na(.x))) %>%
  select(starts_with("SK_L"))

#-------------------------------
# Now we remove the exclusive protein groups from dat_norm since we do not need
# them for the DE test. They are enriched by definition.

#1. Create vector of rownames to remove

exclusive <-
  c(
    rownames(DS_exclusive),
    rownames(DSK_exclusive),
    rownames(K_exclusive),
    rownames(SK_exclusive)
  )

dat_norm_filter <- dat_norm[!(rownames(dat_norm) %in% exclusive), ]


##### Heatmap#################################################
#Generate annotation table
heatanno <- data.frame(sample = colnames(dat.nobatch))
heatanno$Genotype <- str_sub(heatanno$sample, 1, 2)
heatanno$Age <- rep(c(6, 13, 17), 6)
rownames(heatanno) <- colnames(dat.nobatch)
heatanno <- subset (heatanno, select = -sample)
heatanno$Genotype <- factor(heatanno$Genotype, levels = c("WT", "KO"))
heatanno$Age <- factor(heatanno$Age, levels = c(6, 13, 17))

# List with colors for each annotation.
heatcolors <- list(Age = brewer.pal(3, "Set1"),
                   Genotype = brewer.pal(2, "Set2"))
names(heatcolors$Age) <- c("6", "13", "17")
names(heatcolors$Genotype) <- c("WT", "KO")
heatcolors$Genotype <- heatcolors$Genotype[-3]

plotheat <- pheatmap(
  dat.nobatch,
  clustering_distance_cols = "manhattan",
  show_rownames = FALSE,
  cutree_cols = 3,
  color = colorRampPalette(c("navy", "white", "red"))(50),
  na_col = "grey39",
  annotation_col = heatanno,
  annotation_colors = heatcolors
)

png(
  filename = paste0(
    "Heatmap_",
    "log2_abundance_CycL_nobatch_",
    projectname,
    "_500x450px_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  width = 500 ,
  height = 450,
  units = "px"
)
plotheat
dev.off()

# Correlation heatmap###########################
#Compute correlation matrix
cormat <- round(cor(x = dat.nobatch, use = "pairwise.complete.obs"), 2)

#Melt correlation matrix to table which is compatible with ggplot2
melted_cormat <- as.data.frame.table(cormat)


# Get lower triangle of the correlation matrix
get_lower_tri <- function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

#set lower triangle of correlation matrix to NA and save to new matrix
upper_tri <- get_upper_tri(cormat)

# Melt the correlation matrix
melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0.5,
    limit = c(0, 1),
    space = "Lab",
    name = "Pearson\nCorrelation",
    na.value = "white"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      size = 10,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.title = element_text()
  ) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value),
            color = "black",
            size = 3) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.7),
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    barwidth = 7,
    barheight = 1,
    title.position = "top",
    title.hjust = 0.5
  ))

#Save correlation heatmap as png. Adjust image sizes as needed!
ggsave(
  paste0(
    "Correlation_heatmap_",
    "log2_abundance_CycL_noBatch",
    projectname,
    "_200x165mm_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 165,
  width = 200,
  units = "mm",
  dpi = "retina"
)

###### Diferential analysis (LIMMA)


# Make design table to assign sample groups
# if there is only one factor, such as treatment. You can define a vector with
# the treatment group in the same order as samples in the protein table.
cond = as.factor(
  c(
    "WT6",
    "WT13",
    "WT17",
    "KO6",
    "KO13",
    "KO17",
    "WT6",
    "WT13",
    "WT17",
    "KO6",
    "KO13",
    "KO17",
    "WT6",
    "WT13",
    "WT17",
    "KO6",
    "KO13",
    "KO17"
  )
)

# The function model.matrix is used to generate the design matrix. We include batch in model because we see a strong batch effect in PCA/MDS.
# If no batcheffect is present, omit the batch variable.
design = model.matrix( ~ 0 + cond + batch) # 0 means no intercept for the linear model
colnames(design) = gsub("cond", "", colnames(design))


#Make contrast = define comparisons
# you can define one or multiple contrasts here
x <- c(
  "KO6-WT6",
  "KO13-WT13",
  "KO17-WT17",
  "KO17-KO6",
  "WT17-WT6",
  "KO13-KO6",
  "WT13-WT6",
  "KO17-KO13",
  "WT17-WT13"
)
contrast =  makeContrasts(contrasts = x, levels = design)
fit1 <- lmFit(dat_norm, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2, trend = TRUE)



#######################################################
#Save results to plots and tables
for (i in x) {
  limma_result <- topTable(fit3, coef = i, number = Inf)
  limma_result$ProtID <- row.names(limma_result)
  limma_result$GeneID <-
    df_prot[rownames(limma_result), "Gene.Symbol"]
  
  #Add categories "myelin" and "SLI" to limma_results table
  limma_result$category <- ifelse(
    rownames(limma_result) %in% SLI$MasterProtID,
    "SLI",
    ifelse(
      rownames(limma_result) %in% myelin$MasterProtID,
      "myelin",
      "other"
    )
  )
  
  #Reorder dataframe by the new category column
  limma_result <- limma_result %>%
    arrange(factor(category, levels = c("other", "SLI", "myelin")))
  
  # create custom key-value pairs for different cell-types
  # this can be achieved with nested ifelse statements
  keyvals <- ifelse(
    rownames(limma_result) %in% SLI$MasterProtID,
    "darkorange2",
    ifelse(
      rownames(limma_result) %in% myelin$MasterProtID,
      "navyblue",
      "grey"
    )
  )
  keyvals[is.na(keyvals)] <- "grey"
  names(keyvals)[keyvals == 'grey'] <- "other"
  names(keyvals)[keyvals == 'navyblue'] <- "myelin"
  names(keyvals)[keyvals == 'darkorange2'] <- "SLI"
  
  
  ############################
  #Volcano plot
  
  #Generate plot
  EnhancedVolcano(
    limma_result,
    lab = limma_result$GeneID,
    selectLab = limma_result[myelin$MasterProtID, "GeneID"],
    x = 'logFC',
    y = 'P.Value',
    ylim = c(0, 4),
    title = i,
    subtitle = NULL,
    pCutoff = 0.1,
    pCutoffCol = 'adj.P.Val',
    FCcutoff = 0.58,
    max.overlaps = 10,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    endsConnectors = "last",
    legendPosition = "top",
    colAlpha = 0.5,
    colCustom = keyvals,
    caption = bquote( ~ Log[2] ~ "fold change cutoff: 0.58; adj.P.Val cutoff: < 0.10")
  )
  ggsave(
    paste0(
      "Volcano_",
      gsub(" ", "", i),
      "_",
      "CycL_noBatch_adj.P.Val<0.1_customlabel_",
      projectname,
      "_250x180mm_",
      format(Sys.time(), '%Y%m%d_%H%M%S'),
      ".png"
    ),
    height = 250,
    width = 200,
    units = "mm",
    dpi = "retina"
  )
  
  #####################
  #MD-Plot
  limma_plot <- limma_result  %>%
    mutate(
      gene_type = case_when(
        `logFC` > 0.58 & `adj.P.Val` < 0.1 ~ "Up",
        `logFC` < -0.58 &
          `adj.P.Val` < 0.1 ~ "Down",
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
      size = 4,
      force = 1.5,
      nudge_y = 0,
      max.overlaps = 10,
      min.segment.length = 0
    ) +
    scale_fill_manual(values = cols) +  # Modify point colour
    labs(
      title = i,
      x = "Average expression",
      y = "Log2 fold change",
      fill = "Regulation"
    ) +
    theme(
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"),
      legend.title = element_text(face = "bold")
    )
  ggsave(
    paste0(
      "MDplot_",
      gsub(" ", "", i),
      "_",
      "CycL_noBatch_adj.P.Val<0.1",
      projectname,
      "_300x200mm_",
      format(Sys.time(), '%Y%m%d_%H%M%S'),
      ".png"
    ),
    height = 200,
    width = 300,
    units = "mm",
    dpi = "retina"
  )
  
  ########################################
  # Write csv with limma results
  write.csv(
    x = limma_result,
    file = paste0(
      "Limmatrend_",
      i,
      "CycL_noBatch_",
      projectname,
      format(Sys.time(), '%Y%m%d_%H%M%S'),
      ".csv"
    ),
    row.names = FALSE
  )
}

###########################################
## Write output tables for paper
#Candidate tables
candidates6M <- topTable(fit3, coef = "KO6-WT6", number = Inf)
candidates13M <- topTable(fit3, coef = "KO13-WT13", number = Inf)
candidates17M <- topTable(fit3, coef = "KO17-WT17", number = Inf)

candidates6M$Accession <- row.names(candidates6M)
candidates13M$Accession <- row.names(candidates13M)
candidates17M$Accession <- row.names(candidates17M)

#Protein abundance table
prot_abun_filt <- merge(
  x = dat_log_filter,
  y = as.data.frame(dat_norm),
  by = "row.names",
  all = TRUE,
  suffixes = c("_log2", "_norm_log2")
) %>%
  rename(Accession = Row.names)

#Add some additional columns to table from original imported PD result table
prot_abun_filt <- merge(x = prot_abun_filt,
                        y = df_prot[, c(
                          "Accession",
                          "Description",
                          "Gene.Symbol",
                          "Coverage",
                          "Peptides",
                          "PSMs",
                          "Unique.Peptides",
                          "Razor.Peptides",
                          "AA.length",
                          "MW.kDa",
                          "pI",
                          "Mascot.Score",
                          "PEP.Score",
                          "Exp.q_value"
                        )],
                        by = "Accession",
                        all.x = TRUE)

candidates6M <-
  merge(x = candidates6M[, c("Accession", "logFC", "P.Value", "adj.P.Val")],
        y = df_prot[, c("Accession", "Description", "Gene.Symbol")],
        by = "Accession",
        all.x = TRUE)
candidates13M <-
  merge(x = candidates13M[, c("Accession", "logFC", "P.Value", "adj.P.Val")],
        y = df_prot[, c("Accession", "Description", "Gene.Symbol")],
        by = "Accession",
        all.x = TRUE)
candidates17M <-
  merge(x = candidates17M[, c("Accession", "logFC", "P.Value", "adj.P.Val")],
        y = df_prot[, c("Accession", "Description", "Gene.Symbol")],
        by = "Accession",
        all.x = TRUE)

candidates6M$Description <-
  sub("OS.*", "", candidates6M$Description)
candidates13M$Description <-
  sub("OS.*", "", candidates13M$Description)
candidates17M$Description <-
  sub("OS.*", "", candidates17M$Description)

#Add abundance values from normalized matrix for relevant time point
candidates6M <-
  merge(
    x = candidates6M,
    y = select(as.data.frame(dat_norm), contains("6")),
    by.x = "Accession",
    by.y = "row.names"
  )
candidates13M <-
  merge(
    x = candidates13M,
    y = select(as.data.frame(dat_norm), contains("13")),
    by.x = "Accession",
    by.y = "row.names"
  )
candidates17M <-
  merge(
    x = candidates17M,
    y = select(as.data.frame(dat_norm), contains("17")),
    by.x = "Accession",
    by.y = "row.names"
  )

#Reorder column
candidates6M <-
  select(
    candidates6M,
    Gene.Symbol,
    Accession,
    Description,
    logFC,
    P.Value,
    adj.P.Val,
    contains("WT"),
    contains("KO")
  )
candidates13M <-
  select(
    candidates13M,
    Gene.Symbol,
    Accession,
    Description,
    logFC,
    P.Value,
    adj.P.Val,
    contains("WT"),
    contains("KO")
  )
candidates17M <-
  select(
    candidates17M,
    Gene.Symbol,
    Accession,
    Description,
    logFC,
    P.Value,
    adj.P.Val,
    contains("WT"),
    contains("KO")
  )

#Add valid values count per genotype
candidates6M$Observation.count_WT <-
  3 - rowSums(is.na(candidates6M[, 6:8]))
candidates6M$Observation.count_KO <-
  3 - rowSums(is.na(candidates6M[, 9:11]))

candidates13M$Observation.count_WT <-
  3 - rowSums(is.na(candidates13M[, 6:8]))
candidates13M$Observation.count_KO <-
  3 - rowSums(is.na(candidates13M[, 9:11]))

candidates17M$Observation.count_WT <-
  3 - rowSums(is.na(candidates17M[, 6:8]))
candidates17M$Observation.count_KO <-
  3 - rowSums(is.na(candidates17M[, 9:11]))

#Generate filtered tables containing only proteins above adj.P-Value and log2 fold change cutoff
candidates6M.filtered <-
  filter(candidates6M,
         logFC > 0.58 & adj.P.Val < 0.1  | logFC < -0.58 & adj.P.Val < 0.1)
candidates13M.filtered <-
  filter(candidates13M,
         logFC > 0.58 & adj.P.Val < 0.1  | logFC < -0.58 & adj.P.Val < 0.1)
candidates17M.filtered <-
  filter(candidates17M,
         logFC > 0.58 & adj.P.Val < 0.1  | logFC < -0.58 & adj.P.Val < 0.1)

# Write xlsx tables
write.xlsx(
  candidates6M.filtered,
  file = paste0(
    "Candidate_proteins_filtered_",
    "1.5fold_adj.P.Val<0.1_",
    "6M_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".xlsx"
  ),
  keepNA = TRUE
)
write.xlsx(
  candidates13M.filtered,
  file = paste0(
    "Candidate_proteins_filtered_",
    "1.5fold_adj.P.Val<0.1_",
    "13M_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".xlsx"
  ),
  keepNA = TRUE
)
write.xlsx(
  candidates17M.filtered,
  file = paste0(
    "Candidate_proteins_filtered_",
    "1.5fold_adj.P.Val<0.1_",
    "17M_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".xlsx"
  ),
  keepNA = TRUE
)
write.xlsx(
  prot_abun_filt,
  file = paste0(
    "Protein_table_filtered_noNAs_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".xlsx"
  ),
  keepNA = TRUE
)
