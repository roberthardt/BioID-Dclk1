# Cleaning of workspace, console and plot area --------------------------------

# Clear plots
if (!is.null(dev.list()))
  dev.off()
# Clear console
cat("\014")
# Clean workspace
rm(list = ls())

# Define needed package---------------------------------------------------------
package.names <-
  c(
    "BiocManager",
    "tidyverse",
    "tcltk",
    "matrixStats",
    "scales",
    "svglite",
    "ggsci",
    "ggrepel",
    "pheatmap",
    "heatmaply",
    "openxlsx",
    "reshape2",
    "RColorBrewer",
    "imputeLCMD",
    "ggpubr"
  )

# Install/Load packages---------------------------------------------------------

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

BiocManager::install(c('limma', 'EnhancedVolcano', "clusterProfiler",
                       "MsCoreUtils", "pcaMethods", "impute", "org.Mm.eg.db",
                       "NormalyzerDE"))
library("EnhancedVolcano")
library("limma")
library("MsCoreUtils")

# Load data --------------------------------------------------------------------

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

# Set path of working directory
setwd(choose_directory())

# Set output paths for plots and tables
QC_path <- "QC_plots/"

# Select file
psites <-
  read.table(
    file = choose_file(),
    stringsAsFactors = FALSE,
    header = TRUE,
    quote = "",
    comment.char = "",
    sep = "\t"
  )

# Basic data cleaning ----------------------------------------------------------

## Set project name ####
# This name is added to all result tables and plots. 
#Please change for every project/dataset to be analyzed.
projectname <- "03_BioID-DCLK_P-sites"

# Generate sites count table
count_tbl <- data.frame(type = "identified",
                       count = nrow(psites),
                       details = "Number of identified sites in original MQ results txt table")

## Create phosphosite ID column by combining Gene name, Amino acid & Position
psites <- psites %>% unite(col = Site, 
                               c(Leading.proteins, Amino.acid), 
                               remove = FALSE) %>% unite(col = Site,
                                                         c(Site,Position),
                                                         sep = "",
                                                         remove = FALSE)

## Assign unique id column as row names ####
rownames(psites) <- psites$Site

## remove reverse & contaminant entries ####
p_filter <- psites %>%
  filter(Reverse != "+") %>%
  filter(Potential.contaminant != "+")

# Append data to sites count table
count_tbl[nrow(count_tbl)+1,] <- c("identified_wo_rev_cont",
                                  nrow(p_filter),
                                  "Number of identified sites after removal of reverse and contaminants")

## Filter for Class I sites = site probability >0.75 
p_filter <- p_filter %>%
  filter(Localization.prob > 0.75)

# Append data to sites count table
count_tbl[nrow(count_tbl)+1,] <- c("classI",
                                   nrow(p_filter),
                                   "Number of identified class I (loc. prob >0.75) sites")

## Filter for expression columns ####
dat <- select(p_filter, contains(c(
  "Intensity.L.BioID",
  "Intensity.H.BioID"
)))

## Set 0 values to NA ####
dat[dat == 0] <- NA

## Rename columns ####
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

## Log2 transform data ####
dat_log <- dat %>%
  dplyr::mutate(across(.cols = everything(),
                       .fns = log2))

# QC Plots----------------------------------------------
# Generate some basic qc plots (box, density) in ggplot. 
## Reformatting ####
### Therefore convert to long format first####
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


### Reorder factor levels to achieve different sorting in plots####
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

## Make plots####

### Boxplot####

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

# Save boxplot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(QC_path,
    "P-sites_Boxplot_log2_intensities_L+H_400x150_",
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
  paste0(QC_path,
    "P-sites_Boxplot_log2_intensities_L+H_400x150_",
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


### Density plots####
# Calc mean log2 intensities
stat_log2_inten <- dat_log_long %>%
  group_by(label, sampleID) %>%
  summarize(
    mean = mean(log2_intensity, na.rm = TRUE),
    median = median(log2_intensity, na.rm = TRUE)
  )

# Plot with ggplot
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

# Save boxplot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(QC_path,
    "P-sites_Density_log2_intensities_L+H_400x100_",
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
  paste0(QC_path,
    "P-sites_Density_log2_intensities_L+H_400x100_",
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

### MA-Plot####
#Create MA plots against global average
A <- rowMeans(dat_log, na.rm = TRUE)
M <- sweep(dat_log, 1, rowMeans(dat_log, na.rm = TRUE))

pdf(
  file = paste0(QC_path,
    "P-sites_MA-plots_nonorm",
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

# Valid values filtering-------------------------------------------------------------
# Count number of NaN per condition.
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

# Append data to sites count table
count_tbl[nrow(count_tbl)+1,] <- c("2_valid_L_val_anygroup",
                                   nrow(dat_log_filter),
                                   "Sites quantified in light channel for at least 2 times in any sample group")

# Imputation-------------------------------------------------------------------
# Impute missing signals by normal distribution
# First store position info of imputed values in a separate data frame
na_pos <- is.na(dat_log_filter) %>%
  as.data.frame()

# Perform imputation
dat_log_impute <- dat_log_filter %>%
  as.matrix() %>%
  MsCoreUtils::impute_matrix(method = "MinProb") %>%
  as.data.frame()

# QC Plots---------------------------------------------------------------------
# Plot density plots for data after imputation
# Generate some basic qc plots (box, density) in ggplot. 
## Reformatting####
### Therefore convert to long format first
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


## Make plots ####

### Boxplot ####

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
  paste0(QC_path,
    "P-sites_Boxplot_log2_intensities+impute_L+H_400x150_",
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
  paste0(QC_path,
    "P-sites_Boxplot_log2_intensities+impute_L+H_400x150_",
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


### Density plots ####
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

#Save plot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(QC_path,
    "P-sites_Density_log2_intensities+impute_L+H_400x100_",
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
  paste0(QC_path,
    "P-sites_Density_log2_intensities+impute_L+H_400x100_",
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

# Ratio generation and further filtering --------------------------------------
## Generate L/H ratios####
df_ratios <- dat_log_impute[, 1:12] - dat_log_impute[, 13:24]

## Filter for enriched proteins ####
### Count number of L/H ratios >1 per condition.####
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

### Filter for enriched proteins. (at least two log2 ratios > 1)####
df_ratios_filter <- df_ratios %>%
  filter(if_any(contains("upreg_count"),
                ~ . >= 2))
### Generate vector with protein accessions of enriched proteins####
enriched_ID <- row.names(df_ratios_filter)

### Fetch light intensities of enriched proteins BEFORE imputation####
# Use rownames from filtered df_ratios to fetch the corresponding light log2 
# intensities 
dat_log_enriched <- dat_log_filter[enriched_ID,
                                   c(1:12)]

# Append data to sites count table
count_tbl[nrow(count_tbl)+1,] <- c("enriched",
                                   nrow(dat_log_enriched),
                                   "Sites having at least 2 log2 L/H ratios >1 in any sample group")

# QC Plots ---------------------------------------------------------------------
# Plot light log2 intensities in boxplot and density plots of the 
# enriched proteins
## Reformatting#####
### Convert to long format
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

### Reorder factor levels to achieve different sorting in plots
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
## Make plots####
### Boxplot#####

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

# Save boxplot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(QC_path,
    "P-sites_Boxplot_enriched_log2_intensities_L_400x150_",
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
  paste0(QC_path,
    "P-sites_Boxplot_enriched_log2_intensities_L_400x150_",
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


### Density plots####
# Plot with ggplot
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

# Save density plot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(QC_path,
    "P-sites_Density_enriched_log2_intensities_L_180x120_",
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
  paste0(QC_path,
    "P-sites_Density_enriched_log2_intensities_L_180x120_",
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

# Normalization ----------------------------------------------
## Evaluation####
# Check out which normalization method to use best with NormalyzerDE
BiocManager::install("NormalyzerDE", force = TRUE)
library("NormalyzerDE")

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

# Generate input data & write to file.
# Input data has to be de-logarithmized before
input_data <- 2 ^ dat_log_enriched  %>%
  mutate(ProtID = row.names(dat_log_enriched),
         .before = 1) %>%
  write_tsv(file = "inputNormalyzerDE.tsv")

# Run normalyzerDE
NormalyzerDE::normalyzer(
  jobName = "normalizerDE_eval",
  designPath = "designNormalyzerDE.tsv",
  dataPath = "inputNormalyzerDE.tsv",
  outputDir = getwd()
)

# Result: ==> basically every method is better then simple log2 or mean. But RLR,
# cyclicLoess and vsn perform best. Cyclic loess leads to lower Intragroup
# Pearson & Spearman correlations. I opt for RLR.

## Perform normalization ####
dat_norm <-
  performGlobalRLRNormalization(as.matrix(dat_log_enriched), 
                                noLogTransform = TRUE) %>%
  as.data.frame()

# QC Plots--------------------------------------------------------------------
# Plot normalized data in Boxplot, Density plot and PCA

## Reformatting ####
# Convert to long format
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

# Reorder factor levels to achieve different sorting in plots
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

## Make plots####
### Boxplot#####

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
  paste0(QC_path,
    "P-sites_Boxplot_enriched_norm_log2_intensities_L_400x150_",
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
  paste0(QC_path,
    "P-sites_Boxplot_enriched_norm_log2_intensities_L_400x150_",
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

### Density plots ####
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
  paste0(QC_path,
    "P-sites_Density_enriched_norm_log2_intensities_L_180x120_",
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
  paste0(QC_path,
    "P-sites_Density_enriched_norm_log2_intensities_L_180x120_",
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

### PCA####

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
  paste0(QC_path,
    "P-sites_PCA_",
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
  paste0(QC_path,
    "P-sites_PCA_",
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


### Correlation heatmap ####
# Compute correlation matrix
cormat <- round(cor(x = dat_norm, use = "pairwise.complete.obs"), 2)

# Melt correlation matrix to table which is compatible with ggplot2
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

# Set lower triangle of correlation matrix to NA and save to new matrix
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

# Save correlation heatmap as png. Adjust image sizes as needed!
ggsave(
  paste0(QC_path,
    "P-site_Correlation_heatmap_",
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

# Imputation-------------------------------------------------------------------
# We impute missing values of light intensities to get more robust statistics

# Impute missing signals by normal distribution
#First store position info of imputed values in a separate data frame
na_pos_2 <- is.na(dat_norm) %>%
  as.data.frame()

#Perform imputation
dat_norm_impute <- dat_norm %>%
  as.matrix() %>%
  MsCoreUtils::impute_matrix(method = "MinProb") %>%
  as.data.frame()

# QC Plots---------------------------------------------------------------------
#Plot data after imputation
#Generate some basic qc plots (box, density, pca) in ggplot. 

##Reformatting####
#Therefore convert to long format first
dat_norm_impute_long <- dat_norm_impute %>%
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

#Add ID column
dat_norm_impute_long <- mutate(.data = dat_norm_impute_long,
                              id = row.names(dat_norm_impute_long))

#Pivot also the NA-Positions to the long format
na_pos_long_2 <- na_pos_2 %>%
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
  mutate(replicate = factor(replicate, levels = unique(replicate))) %>%
  mutate(sampleID = factor(sampleID, levels = unique(sampleID))) %>%
  mutate(imputed = factor(imputed, levels = unique(imputed)))

#Add ID column
na_pos_long_2 <- mutate(.data = na_pos_long_2,
                      id = row.names(na_pos_long_2))

#Merge expression data frame and na_position frame
dat_norm_impute_merge <- merge(x = dat_norm_impute_long,
                              y = na_pos_long_2[, c("id", "imputed")],
                              by = "id")

#reorder factor levels to achieve different sorting in plots
levels(dat_norm_impute_merge$sample) <- c("DSK", "DS", "SK", "K")
levels(dat_norm_impute_merge$replicate) <- c("01", "02", "03")
levels(dat_norm_impute_merge$sampleID) <-
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


## Make plots####

### Boxplot####

# Calc mean log2 intensities
stat_log2_inten <- dat_norm_impute_merge %>%
  summarize(
    mean = mean(log2_intensity, na.rm = TRUE),
    median = median(log2_intensity, na.rm = TRUE)
  )

# Plot with ggplot
ggplot(dat_norm_impute_merge,
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
  paste0(QC_path,
    "P-sites_Boxplot_enriched_norm_log2_intensities+impute_L_400x150_",
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
  paste0(QC_path,
    "P-sites_Boxplot_enriched_norm_log2_intensities+impute_L_400x150_",
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


### Density plots####
# Calc mean log2 intensities
stat_log2_inten <- dat_norm_impute_merge %>%
  group_by(sampleID) %>%
  summarize(
    mean = mean(log2_intensity, na.rm = TRUE),
    median = median(log2_intensity, na.rm = TRUE)
  )

#Plot with ggplot
ggplot(dat_norm_impute_merge, aes(x = log2_intensity,
                                 fill = imputed)) +
  geom_density(alpha = 0.75)  +
  facet_wrap(~ sampleID,
             ncol = 3) +
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

#Save plot as png and svg. Adjust image sizes as needed!
ggsave(
  paste0(QC_path,
    "P-sites_Density_enriched_norm_log2_intensities+impute_L_200x200_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  height = 200,
  width = 200,
  units = "mm",
  dpi = "retina"
)

ggsave(
  paste0(QC_path,
    "P-sites_Density_enriched_norm_log2_intensities+impute_L_200x200_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".svg"
  ),
  height = 200,
  width = 200,
  units = "mm",
  dpi = "retina"
)

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
    into = c("Sample_type", "Label", "Replicate"),
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
  paste0(QC_path,
    "P-sites_PCA_",
    "enriched_norm_log2_intensities+impute_L",
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
    "enriched_norm_log2_intensities+impute_L",
    "_200x165mm_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".svg"
  ),
  height = 165,
  width = 200,
  units = "mm",
  dpi = "retina"
)

### Heatmap####
#Generate annotation table
heatanno <- data.frame(sample = colnames(dat_norm_impute))
heatanno$Construct <- str_extract(heatanno$sample, "[^_]+")
rownames(heatanno) <- colnames(dat_norm_impute)
heatanno <- subset (heatanno, select = -sample)
heatanno$Construct <- factor(heatanno$Construct, levels = c("DSK", "DS", "SK", "K"))


# List with colors for each annotation.
heatcolors <- list(Construct = brewer.pal(4, "Dark2"))
names(heatcolors$Construct) <- c("DSK", "DS", "SK", "K")

plotheat <- pheatmap(
  dat_norm_impute,
  clustering_distance_cols = "manhattan",
  show_rownames = FALSE,
  cutree_cols = 4,
  color = colorRampPalette(c("navy", "white", "red"))(50),
  na_col = "grey39",
  annotation_col = heatanno,
  annotation_colors = heatcolors
)

png(
  filename = paste0(QC_path,
    "P-sites_Heatmap_",
    "enriched_norm_log2_intensities+impute_L",
    "_1000x700px_",
    projectname,
    "_",
    format(Sys.time(), '%Y%m%d_%H%M%S'),
    ".png"
  ),
  width = 1000 ,
  height = 700,
  units = "px"
)
plotheat
dev.off()

# Differential analysis (LIMMA)-------------------------------------------------
## Calculations####
# Make design table to assign sample groups
# if there is only one factor, such as treatment. You can define a vector with
# the treatment group in the same order as samples in the protein table.
cond = as.factor(
  c(
    "DS",
    "DS",
    "DS",
    "DSK",
    "DSK",
    "DSK",
    "K",
    "K",
    "K",
    "SK",
    "SK",
    "SK"
  )
)

# The function model.matrix is used to generate the design matrix. 
# We should include batch in model if we see a strong batch effect in PCA/MDS.
# If no batcheffect is present, omit the batch variable.
design = model.matrix( ~ 0 + cond) # 0 means no intercept for the linear model
colnames(design) = gsub("cond", "", colnames(design))

# Make contrast = define comparisons
# you can define one or multiple contrasts here
x <- c(
  "DSK-DS",
  "DSK-SK",
  "DSK-K",
  "DS-SK",
  "DS-K",
  "SK-K"
)
contrast =  makeContrasts(contrasts = x, levels = design)
fit1 <- lmFit(as.matrix(dat_norm_impute), design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2, trend = TRUE)


## Save results to plots and tables####

DE_path <- "DE_analysis/" # Set path for output data

for (i in x) {
  limma_result <- topTable(fit3, coef = i, number = Inf)
  limma_result$SiteID <- row.names(limma_result)
  limma_result$GeneID <- psites[row.names(limma_result), "Gene.names"]
  limma_result$AA <- psites[rownames(limma_result), "Amino.acid"]
  limma_result$Position <- psites[rownames(limma_result), "Position"]
  limma_result <- limma_result %>% 
    unite(col = Site, c(GeneID, AA)) %>%
    unite(col = Site,c(Site, Position), sep = "")
  
  ### Volcano plot####
  
  # Generate plot
  EnhancedVolcano(
    limma_result,
    lab = limma_result$Site,
    x = 'logFC',
    y = 'P.Value',
    title = i,
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
      "P-sites_Volcano_",
      gsub(" ", "", i),
      "_",
      "Limma_norm_imputed_adj.P.Val<0.05_",
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
  
  ## MD-Plot#####
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
      aes(label = Site),
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
      DE_path,
      "P-sites_MDplot_",
      gsub(" ", "", i),
      "_",
      "Limma_norm_imputed_adj.P.Val<0.05_",
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
  
  ## Write csv with limma results####
  write.csv(
    x = limma_result,
    file = paste0(
      DE_path,
      "P-sites_Limma_results_table_",
      i,
      "_norm_imputed_",
      projectname,
      format(Sys.time(), '%Y%m%d_%H%M%S'),
      ".csv"
    ),
    row.names = FALSE
  )
}

# Write output tables for paper------------------------------------------------
## Candidate tables####
candidates_DSK_DS <- topTable(fit3, coef = "DSK-DS", number = Inf) %>%
  dplyr::select(c("logFC", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(`logFC` > 0.58 & `adj.P.Val` < 0.05 |
                  `logFC` < -0.58 & `adj.P.Val` < 0.05 ) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_norm %>%
              as.data.frame() %>%
              select(starts_with(c("DSK_","DS_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_norm_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>%
  left_join(x = ., 
            y = psites %>%
              select(c("Site", "Gene.names", "Amino.acid", "Position", "Fasta.headers")) %>%
              rename(SiteID = Site))

candidates_DSK_SK <- topTable(fit3, coef = "DSK-SK", number = Inf) %>%
  dplyr::select(c("logFC", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(`logFC` > 0.58 & `adj.P.Val` < 0.05 |
                  `logFC` < -0.58 & `adj.P.Val` < 0.05 ) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_norm %>%
              as.data.frame() %>%
              select(starts_with(c("DSK_","SK_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_norm_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>%
  left_join(x = ., 
            y = psites %>%
              select(c("Site", "Gene.names", "Amino.acid", "Position", "Fasta.headers")) %>%
              rename(SiteID = Site))

candidates_DSK_K <- topTable(fit3, coef = "DSK-K", number = Inf) %>%
  dplyr::select(c("logFC", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(`logFC` > 0.58 & `adj.P.Val` < 0.05 |
                  `logFC` < -0.58 & `adj.P.Val` < 0.05 ) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_norm %>%
              as.data.frame() %>%
              select(starts_with(c("DSK_","K_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_norm_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>%
  left_join(x = ., 
            y = psites %>%
              select(c("Site", "Gene.names", "Amino.acid", "Position", "Fasta.headers")) %>%
              rename(SiteID = Site))

candidates_DS_SK <- topTable(fit3, coef = "DS-SK", number = Inf) %>%
  dplyr::select(c("logFC", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(`logFC` > 0.58 & `adj.P.Val` < 0.05 |
                  `logFC` < -0.58 & `adj.P.Val` < 0.05 ) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_norm %>%
              as.data.frame() %>%
              select(starts_with(c("DS_","SK_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_norm_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>%
  left_join(x = ., 
            y = psites %>%
              select(c("Site", "Gene.names", "Amino.acid", "Position", "Fasta.headers")) %>%
              rename(SiteID = Site))

candidates_DS_K <- topTable(fit3, coef = "DS-K", number = Inf) %>%
  dplyr::select(c("logFC", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(`logFC` > 0.58 & `adj.P.Val` < 0.05 |
                  `logFC` < -0.58 & `adj.P.Val` < 0.05 ) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_norm %>%
              as.data.frame() %>%
              select(starts_with(c("DS_","K_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_norm_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>%
  left_join(x = ., 
            y = psites %>%
              select(c("Site", "Gene.names", "Amino.acid", "Position", "Fasta.headers")) %>%
              rename(SiteID = Site))

candidates_SK_K <- topTable(fit3, coef = "SK-K", number = Inf) %>%
  dplyr::select(c("logFC", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(`logFC` > 0.58 & `adj.P.Val` < 0.05 |
                  `logFC` < -0.58 & `adj.P.Val` < 0.05 ) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_norm %>%
              as.data.frame() %>%
              select(starts_with(c("SK_","K_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_norm_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>%
  left_join(x = ., 
            y = psites %>%
              select(c("Site", "Gene.names", "Amino.acid", "Position", "Fasta.headers")) %>%
              rename(SiteID = Site))

## Abundance table####
#Generate p-site table with abundance values from normalized matrix and 
#  norm+imputed matrix
# First add the normalized abundances
psite_table <-
  merge(
    x = psites %>%
      dplyr::select(c("Site",
                      "Leading.proteins",
                      "Gene.names",
                      "Amino.acid",
                      "Position",
                      "Fasta.headers",
                      "Localization.prob",
                      "Score",
                      "PEP")),
    y = dat_norm %>%
      as.data.frame() %>%
      rename_with(.cols = everything(), 
                  function(x){paste0("log2_norm_intensity_", x)}) %>%
      rownames_to_column(var = "Site"),
    by = "Site"
  )

#Then also add the imputed abundances
psite_table <-
  left_join(
    x = psite_table,
    y = dat_norm_impute %>%
      as.data.frame() %>%
      rename_with(.cols = everything(), 
                  function(x){paste0("log2_norm_imputed_intensity_", x)}) %>%
      rownames_to_column(var = "Site"),
    by = "Site"
  )

#Add indicator columns marking imputed values
psite_table <-
  left_join(
    x = psite_table,
    y = na_pos_2 %>% 
      rename_with(.cols = everything(), 
                  function(x){paste0("is_imputed_", x)}) %>%
      rownames_to_column(var = "Site"),
    by = "Site"
  )

## Combined output workbook####
# Write xlsx workbook which combines prot_table and the candidate tables
wbook <- createWorkbook("P-sites_Results")
addWorksheet(wbook, "Enriched_P-sites")
addWorksheet(wbook, "Candidates_DSKvsDS")
addWorksheet(wbook, "Candidates_DSKvsSK")
addWorksheet(wbook, "Candidates_DSKvsK")
addWorksheet(wbook, "Candidates_DSvsSK")
addWorksheet(wbook, "Candidates_DSvsK")
addWorksheet(wbook, "Candidates_SKvsK")
writeData(wbook, sheet = 1, psite_table)
writeData(wbook, sheet = 2, candidates_DSK_DS)
writeData(wbook, sheet = 3, candidates_DSK_SK)
writeData(wbook, sheet = 4, candidates_DSK_K)
writeData(wbook, sheet = 5, candidates_DS_SK)
writeData(wbook, sheet = 6, candidates_DS_K)
writeData(wbook, sheet = 7, candidates_SK_K)
saveWorkbook(wbook,
             paste0(DE_path, "Enriched_P-sites_and_DE_results_",
                    projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".xlsx"),
             overwrite = TRUE)


# Dclk1 detailed psite analysis ####----------------------------------------
## Filter for Dclk1 psites #####
Dclk1_sites <- psite_table %>% 
  filter(grepl("Q9JLM8",Leading.proteins))

## Create new site label ####
Dclk1_sites <- Dclk1_sites %>% 
  unite(col = Site2, c(Amino.acid, Position),
        remove = FALSE)

## Select only site annotation and norm_intensity columns ####
Dclk1_sites_ma <- select(Dclk1_sites, 
                         c(Position, 
                           Site2, 
                           contains("norm_intensity"))) %>% 
  dplyr:::arrange(Position)


## Set rownames needed for heatmap ####
row.names(Dclk1_sites_ma) <- Dclk1_sites_ma$Site2

## Remove unnecessary columns ####
Dclk1_sites_ma <- Dclk1_sites_ma %>% 
  select(contains("norm_intensity"))

### Set colnames needed for pheatmap col_labels ####
colnames(Dclk1_sites_ma) <- c("DS_L_01",
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
                            "SK_L_03")

### Arrange columns in specific order ####
Dclk1_sites_ma <- Dclk1_sites_ma %>%
  select("DSK_L_01",
         "DSK_L_02",
         "DSK_L_03",
         "DS_L_01",
         "DS_L_02",
         "DS_L_03",
         "SK_L_01",
         "SK_L_02",
         "SK_L_03",
         "K_L_01",
         "K_L_02",
         "K_L_03")

### Remove DS samples ####
Dclk1_sites_ma_2 <- Dclk1_sites_ma %>%
  select("DSK_L_01",
         "DSK_L_02",
         "DSK_L_03",
         "SK_L_01",
         "SK_L_02",
         "SK_L_03",
         "K_L_01",
         "K_L_02",
         "K_L_03")

### Generate Heatmap####
#### Generate annotation table####
heatanno <- data.frame(sample = colnames(Dclk1_sites_ma))
heatanno$Construct <- rep(c("DSK", "DS", "SK", "K"),
                          each = 3)
rownames(heatanno) <- colnames(Dclk1_sites_ma)
heatanno <- subset (heatanno, select = -sample)
heatanno$Construct <- factor(heatanno$Construct, levels = c("DSK", "DS", "SK", "K"))


#### List with colors for each annotation####
heatcolors <- list(Construct = brewer.pal(4, "Dark2"))
names(heatcolors$Construct) <- c("DSK", "DS", "SK", "K")

#### List with colors without group DS####
heatcolors2 <- list(Construct = heatcolors$Construct[-2])

#### Plot & Save heatmap - all#####
plotheat <- pheatmap(
  Dclk1_sites_ma,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cutree_cols = 4,
  color = colorRampPalette(c("navy", "white", "red"))(50),
  na_col = "grey39",
  scale = "row",
  annotation_col = heatanno,
  annotation_colors = heatcolors,
  fontsize = 16,
  angle_col = 45,
  cellwidth = 30
)

png(
  filename = paste0(QC_path,
                    "Dclk1_Sitemap_",
                    "enriched_norm_log2_intensities_L_z-scaled",
                    "_10x8in_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".png"
  ),
  width = 10 ,
  height = 8,
  units = "in",
  res = 300
)
plotheat
dev.off()

 svg(
  filename = paste0(QC_path,
                    "Dclk1_Sitemap_",
                    "enriched_norm_log2_intensities_z-scaled",
                    "_8x6in_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".svg"
  ),
  width = 10 ,
  height = 8
)
plotheat
dev.off()


#### Plot & Save heatmap - w/o DS#####
plotheat <- pheatmap(
  Dclk1_sites_ma_2,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cutree_cols = 4,
  color = colorRampPalette(c("navy", "white", "red"))(50),
  na_col = "grey39",
  scale = "row",
  annotation_col = heatanno,
  annotation_colors = heatcolors2,
  fontsize = 16,
  angle_col = 45,
  cellwidth = 30
)

png(
  filename = paste0(QC_path,
                    "Dclk1_Sitemap_",
                    "enriched_norm_log2_intensities_L_z-scaled_woDS",
                    "_10x8in_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".png"
  ),
  width = 10 ,
  height = 8,
  units = "in",
  res = 300
)
plotheat
dev.off()

svg(
  filename = paste0(QC_path,
                    "Dclk1_Sitemap_",
                    "enriched_norm_log2_intensities_z-scaled_woDS",
                    "_8x6in_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".svg"
  ),
  width = 10 ,
  height = 8
)
plotheat
dev.off()


# Regulated psites heatmaps####-----------------------------------------------
## Filter prot_table for all regulated psites ####
psite_regulated <- psite_table %>% 
  filter(Site %in% 
           c(candidates_DSK_DS$SiteID,
             candidates_DSK_K$SiteID,
             candidates_DSK_SK$SiteID,
             candidates_DS_K$SiteID,
             candidates_DS_SK$SiteID,
             candidates_SK_K$SiteID))

# Append data to sites count table
count_tbl[nrow(count_tbl)+1,] <- c("regulated",
                                   nrow(psite_regulated),
                                   "Sites showing significant regulation in any pairwise group comparison")

## Generate count plot for regulated sites ####
### Convert to long format ####
psite_reg_long <- psite_regulated %>%
  pivot_longer(cols = contains("log2_norm_intensity"),
               names_to = c("Sample","Rep"),
               values_to ="log2_norm_intensity",
               names_sep = "_L_") %>%
  separate(
    col = "Sample",
    into = c("Abundance_type", "Sample_type"),
    sep = "log2_norm_intensity_"
  )

### Generate summary table ####
psite_reg_count <- psite_reg_long %>%
  dplyr::select(Sample_type, Rep, log2_norm_intensity) %>%
  dplyr::group_by(Sample_type, Rep) %>%
  dplyr::summarise(Count = sum(!is.na(log2_norm_intensity))) %>%
  dplyr::mutate(Sample_type = factor(Sample_type, 
                                     levels = c("DSK","DS","SK", "K")))

### Plot barplot ####
ggpubr::ggbarplot(psite_reg_count, 
                  x = "Sample_type", y = "Count", 
                  add = c("mean_sd", "jitter"),
                  fill = "Sample_type",
                  ncol = 1,
                  position = position_dodge(0.8),
                  size = 1,
                  alpha = 0.8,
                  palette = "Dark2",
                  legend = "none",
                  x.text.angle = 45,
                  font.x = c(18, "bold"),
                  font.y = c(18, "bold"),
                  font.main = c(18, "bold"),
                  font.tickslab = c(16, "bold"),
                  title = "Class I phosphosite count",
                  subtitle = "only regulated sites found >=2 in one condition ",
                  xlab = "Construct",
                  ylab = "Count")

ggsave(
  paste0(QC_path,
         "Barplot_validvalcount_",
         "enriched_norm_log2_intensities_L_only_regulated",
         "_120x150px_",
         projectname,
         "_",
         format(Sys.time(), '%Y%m%d_%H%M%S'),
         ".png"
  ),
  width = 120,
  height = 150,
  units = "mm",
  dpi = "retina"
)

## Generate UpSet plot for regulated psites ####
### Generate logical matrix ####
psite_reg_upset <- psite_regulated %>%
  dplyr::select(contains("norm_intensity"))
psite_reg_upset <- !is.na(psite_reg_upset)
psite_reg_upset <- psite_reg_upset %>%
  as.data.frame() %>%
  mutate_all(as.integer)

###Plot UpSet ####
png(
  filename = paste0(QC_path,
                    "Upset_",
                    "enriched_norm_log2_intensities_L_only_regulated",
                    "_1400x800px_",
                    projectname, 
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".png"
  ),
  width = 1400 ,
  height = 800,
  units = "px"
)
upset(psite_reg_upset, 
      nsets = 12,
      text.scale = 2)
dev.off()

## Create heatmap for all regulated psites ####

### Extract only intensities####

#### w/o imputation
psite_reg_ma <-psite_regulated %>%
  select(contains("norm_intensity"))

#### with imputation
psite_reg_ma_impute <-psite_regulated %>%
  select(contains("norm_imputed_intensity"))

### Set rownames for matrix####
rownames(psite_reg_ma) <- psite_regulated$Site
rownames(psite_reg_ma_impute) <- psite_regulated$Site

### Set colnames needed for pheatmap col_labels ####
colnames(psite_reg_ma) <- c("DS_L_01",
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
                               "SK_L_03")

colnames(psite_reg_ma_impute) <- c("DS_L_01",
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
                            "SK_L_03")

### Generate Heatmap####
####Generate annotation table####
heatanno <- data.frame(sample = colnames(psite_reg_ma))
heatanno$Construct <- rep(c("DS", "DSK", "K", "SK"),
                          each = 3)
rownames(heatanno) <- colnames(psite_reg_ma)
heatanno <- subset (heatanno, select = -sample)
heatanno$Construct <- factor(heatanno$Construct, levels = c("DSK", "DS", "SK", "K"))

#### List with colors for each annotation####
heatcolors <- list(Construct = brewer.pal(4, "Dark2"))
names(heatcolors$Construct) <- c("DSK", "DS", "SK", "K")


#### List with colors without group DS####
heatcolors2 <- list(Construct = heatcolors$Construct[-2])

#### Map with imputations - all #####
plotheat_impute <- pheatmap(
  psite_reg_ma_impute,
  clustering_distance_cols = "manhattan",
  show_rownames = F,
  cutree_cols = 4,
  color = colorRampPalette(c("navy", "white", "red"))(50),
  na_col = "grey39",
  scale = "row",
  annotation_col = heatanno,
  annotation_colors = heatcolors,
  fontsize = 16,
  angle_col = 45 
)

#### Extract ordering of all heatmap ####
order_all <- plotheat_impute$tree_row[["order"]]

#### Save map ####
png(
  filename = paste0(QC_path,
                    "Heatmap_",
                    "enriched_norm_log2_intensities+impute_L_only_regulated_z-scaled",
                    "_800x600px_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".png"
  ),
  width = 800 ,
  height = 600,
  units = "px"
)
plotheat_impute
dev.off()

svg(
  filename = paste0(QC_path,
                    "Heatmap_",
                    "enriched_norm_log2_intensities+impute_L_only_regulated_z-scaled",
                    "_8x6in_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".svg"
  ),
  width = 8 ,
  height = 6
)
plotheat_impute
dev.off()

#### Map with imputations - w/o DS #####
plotheat_impute <- pheatmap(
  psite_reg_ma_impute[,4:12],
  clustering_distance_cols = "manhattan",
  show_rownames = F,
  cutree_cols = 3,
  color = colorRampPalette(c("navy", "white", "red"))(50),
  na_col = "grey39",
  scale = "row",
  annotation_col = heatanno,
  annotation_colors = heatcolors2,
  fontsize = 16,
  angle_col = 45 
)

#### Extract ordering of w/o DS heatmap ####
order_woDS <- plotheat_impute$tree_row[["order"]]

#### Save map ####
png(
  filename = paste0(QC_path,
                    "Heatmap_",
                    "enriched_norm_log2_intensities+impute_L_only_regulated_z-scaled_woDS",
                    "_800x600px_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".png"
  ),
  width = 800 ,
  height = 600,
  units = "px"
)
plotheat_impute
dev.off()

svg(
  filename = paste0(QC_path,
                    "Heatmap_",
                    "enriched_norm_log2_intensities+impute_L_only_regulated_z-scaled_woDS",
                    "_8x6in_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".svg"
  ),
  width = 8 ,
  height = 6
)
plotheat_impute
dev.off()


#### Map without imputations - all #####
plotheat <- pheatmap(
  psite_reg_ma[order_all,], # use ordering from clustered all heatmap
  cluster_rows = FALSE,
  clustering_distance_cols = "manhattan",
  show_rownames = F,
  cutree_cols = 4,
  color = colorRampPalette(c("navy", "white", "red"))(50),
  na_col = "grey39",
  scale = "row",
  annotation_col = heatanno,
  annotation_colors = heatcolors,
  fontsize = 16,
  angle_col = 90
)

png(
  filename = paste0(QC_path,
                    "Heatmap_",
                    "enriched_norm_log2_intensities_L_only_regulated_z-scaled",
                    "_800x600px_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".png"
  ),
  width = 800 ,
  height = 600,
  units = "px"
)
plotheat
dev.off()

svg(
  filename = paste0(QC_path,
                    "Heatmap_",
                    "enriched_norm_log2_intensities_L_only_regulated_z-scaled",
                    "_8x6in_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".svg"
  ),
  width = 8 ,
  height = 6
)
plotheat
dev.off()

#### Map without imputations - w/o DS#####
plotheat <- pheatmap(
  psite_reg_ma[order_woDS,4:12],
  cluster_rows = FALSE,
  clustering_distance_cols = "manhattan",
  show_rownames = F,
  cutree_cols = 3,
  color = colorRampPalette(c("navy", "white", "red"))(50),
  na_col = "grey39",
  scale = "row",
  annotation_col = heatanno,
  annotation_colors = heatcolors2,
  fontsize = 16,
  angle_col = 90
)

png(
  filename = paste0(QC_path,
                    "Heatmap_",
                    "enriched_norm_log2_intensities_L_only_regulated_z-scaled_woDS",
                    "_800x600px_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".png"
  ),
  width = 800 ,
  height = 600,
  units = "px"
)
plotheat
dev.off()

svg(
  filename = paste0(QC_path,
                    "Heatmap_",
                    "enriched_norm_log2_intensities_L_only_regulated_z-scaled_woDS",
                    "_8x6in_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".svg"
  ),
  width = 8 ,
  height = 6
)
plotheat
dev.off()


# Write site id count table to Excel ####
wbook <- createWorkbook("P-sites_ID-Count_table")
addWorksheet(wbook, "ID-counts")
writeData(wbook, sheet = 1, count_tbl)
saveWorkbook(wbook,
             paste0(DE_path, "P-sites_ID-count_table_",
                    projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".xlsx"),
             overwrite = TRUE)

# Annotation-enrichment analysis-------------------------------------------------------
# Perform GO-overrepresentation and GSEA analyses with the clusterProfiler 
# package

## Load packages#####
library(clusterProfiler) #load clusterProfiler package
library(org.Mm.eg.db) # load annotation database for Mus musculus
library(DOSE) # load DOSE package for parsing results
library(enrichplot) # needed for GO plots

## Create output subfolder
if (dir.exists("Annotation_enrichment")) {
  
  cat("The folder already exists")
  
} else {
  
  dir.create("Annotation_enrichment")
  
}

## Set relative output path####
annot_path <- "Annotation_enrichment/"

## Optional - Load previous results ####
# If the upper part of teh script, e.g. the DE enrichment analyis,
# was performed before we can also reload the data from the excel
# file here

# Load workbook
wbook <- openxlsx::loadWorkbook(file = choose_file())

#Re-extract data frames
psite_table <- openxlsx::read.xlsx(xlsxFile = wbook, sheet = 1)
candidates_DSK_DS <- openxlsx::read.xlsx(xlsxFile = wbook, sheet = 2)
candidates_DSK_SK <- openxlsx::read.xlsx(xlsxFile = wbook, sheet = 3)
candidates_DSK_K <- openxlsx::read.xlsx(xlsxFile = wbook, sheet = 4)
candidates_DS_SK <- openxlsx::read.xlsx(xlsxFile = wbook, sheet = 5)
candidates_DS_K <- openxlsx::read.xlsx(xlsxFile = wbook, sheet = 6)
candidates_SK_K <- openxlsx::read.xlsx(xlsxFile = wbook, sheet = 7)

## GO-terms####
# First we generate the background lists (parameter: "universe") of all protein ids 
# followed by the individual set of differential enriched proteins for every 
# comparison (parameter: "gene")

### Over-representation analysis ####
#### Generate background list####
background_prot_ids <- psite_table$Gene.names %>%
  str_split(";", simplify = TRUE) %>%
  as.data.frame()  %>% 
  dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
background_prot_ids <- background_prot_ids$V1
background_prot_ids <- na.omit(background_prot_ids)
background_prot_ids <- unique(background_prot_ids)

#### Generate DE-protein lists####
y <- c("candidates_DSK_DS", "candidates_DSK_SK", "candidates_DSK_K",
       "candidates_DS_SK", "candidates_DS_K", "candidates_SK_K")

for (i in y) {
  up_ids <- filter(get(i),
               logFC > 0.58 & adj.P.Val < 0.05 ) %>%
    dplyr::select(c("Gene.names"))
  up_ids <- str_split(up_ids$Gene.names, ";", simplify = TRUE) %>%
    as.data.frame() %>% 
    dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
  up_ids <- up_ids$V1
  up_ids <- na.omit(up_ids)
  up_ids <- unique(up_ids)
  
  down_ids <- filter(get(i),
               logFC < -0.58 & adj.P.Val < 0.05 ) %>%
    dplyr::select(c("Gene.names"))
  down_ids <- str_split(down_ids$Gene.names, ";", simplify = TRUE) %>%
    as.data.frame() %>% 
    dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
  down_ids <- down_ids$V1
  down_ids <- na.omit(down_ids)
  down_ids <- unique(down_ids)

#### Perform over-representation analyses####
  ego_up_BP <- enrichGO(
    gene          = up_ids,
    universe      = background_prot_ids,
    keyType = "SYMBOL",
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.2,
    readable      = FALSE
  )

  ego_up_CC <- enrichGO(
    gene          = up_ids,
    universe      = background_prot_ids,
    keyType = "SYMBOL",
    OrgDb         = org.Mm.eg.db,
    ont           = "CC",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.2,
    readable      = FALSE
  )
  ego_up_MF <- enrichGO(
    gene          = up_ids,
    universe      = background_prot_ids,
    keyType = "SYMBOL",
    OrgDb         = org.Mm.eg.db,
    ont           = "MF",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.2,
    readable      = FALSE
  )
  
  ego_down_BP <- enrichGO(
    gene          = down_ids,
    universe      = background_prot_ids,
    keyType = "SYMBOL",
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.2,
    readable      = FALSE
  )
  
  ego_down_CC <- enrichGO(
    gene          = down_ids,
    universe      = background_prot_ids,
    keyType = "SYMBOL",
    OrgDb         = org.Mm.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.2,
    readable      = FALSE
  )
  
  ego_down_MF <- enrichGO(
    gene          = down_ids,
    universe      = background_prot_ids,
    keyType = "SYMBOL",
    OrgDb         = org.Mm.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.2,
    readable      = FALSE
  )

  
  #### Calculate fold enrichment and add it as an extra column####
  ego_up_BP@result <- ego_up_BP@result %>% 
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  
  ego_up_CC@result <- ego_up_CC@result %>% 
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  
  ego_up_MF@result <- ego_up_MF@result %>% 
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  
  ego_down_BP@result <- ego_down_BP@result %>% 
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  
  ego_down_CC@result <- ego_down_CC@result %>% 
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  
  ego_down_MF@result <- ego_down_MF@result %>% 
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  
  
  #### Simplify terms and plot####
  # GO-Terms are simplified by removing redundant terms (semantic similarity)
  # Afterwards plots are generated if enriched terms are found
  
  if(nrow(base::as.data.frame(ego_up_BP@result)) > 0)
  {
    ego_up_BP <- enrichplot::pairwise_termsim(ego_up_BP)
    ego_up_BP_simple <- clusterProfiler::simplify(ego_up_BP, 
                                                  cutoff = 0.7, 
                                                  by = "p.adjust", 
                                                  select_fun = min)
    # Bubble plot
    dotplot(ego_up_BP_simple, 
            x = "FoldEnrichment", 
            title = "Up - Biological Process", 
            showCategory = 30)
    ggsave(
      paste0(
        annot_path,
        "GO_Bubble_up_BP_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top30_",
        projectname,
        "_200x250mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 250,
      width = 200,
      units = "mm",
      dpi = "retina"
    )
  } else
  {
    ego_up_BP_simple <- NULL
  }
  
  if(nrow(base::as.data.frame(ego_up_CC)) > 0)
  {
    ego_up_CC <- pairwise_termsim(ego_up_CC)
    ego_up_CC_simple <- clusterProfiler::simplify(ego_up_CC, 
                                                  cutoff = 0.7, 
                                                  by = "p.adjust", 
                                                  select_fun = min)
    # Bubble plot
    dotplot(ego_up_CC_simple, 
            x = "FoldEnrichment", 
            title = "Up - Cellular Component", 
            showCategory = 30)
    ggsave(
      paste0(
        annot_path,
        "GO_Bubble_up_CC_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top30_",
        projectname,
        "_200x250mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 250,
      width = 200,
      units = "mm",
      dpi = "retina"
    )
  } else
  {
    ego_up_CC_simple <- NULL
  }
  
  if(nrow(base::as.data.frame(ego_up_MF)) > 0)
  {
    ego_up_MF <- pairwise_termsim(ego_up_MF)
    ego_up_MF_simple <- clusterProfiler::simplify(ego_up_MF, 
                                                  cutoff = 0.7, 
                                                  by = "p.adjust", 
                                                  select_fun = min)
    # Bubble plot
    dotplot(ego_up_MF_simple, 
            x = "FoldEnrichment", 
            title = "Up - Molecular Function", 
            showCategory = 30)
    ggsave(
      paste0(
        annot_path,
        "GO_Bubble_up_MF_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top30_",
        projectname,
        "_200x250mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 250,
      width = 200,
      units = "mm",
      dpi = "retina"
    )
  } else
  {
    ego_up_MF_simple <- NULL
  }
  
  if(nrow(base::as.data.frame(ego_down_BP)) > 0)
  {
    ego_down_BP <- pairwise_termsim(ego_down_BP)
    ego_down_BP_simple <- clusterProfiler::simplify(ego_down_BP, 
                                                  cutoff = 0.7, 
                                                  by = "p.adjust", 
                                                  select_fun = min)
    # Bubble plot
    dotplot(ego_down_BP_simple, 
            x = "FoldEnrichment", 
            title = "Down - Biological Process", 
            showCategory = 30)
    ggsave(
      paste0(
        annot_path,
        "GO_Bubble_down_BP_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top30_",
        projectname,
        "_200x250mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 250,
      width = 200,
      units = "mm",
      dpi = "retina"
    )
  } else
  {
    ego_down_BP_simple <- NULL
  }
  
  if(nrow(base::as.data.frame(ego_down_CC)) > 0)
  {
    ego_down_CC <- pairwise_termsim(ego_down_CC)
    ego_down_CC_simple <- clusterProfiler::simplify(ego_down_CC, 
                                                    cutoff = 0.7, 
                                                    by = "p.adjust", 
                                                    select_fun = min)
    # Bubble plot
    dotplot(ego_down_CC_simple, 
            x = "FoldEnrichment", 
            title = "Down - Cellular Component", 
            showCategory = 30)
    ggsave(
      paste0(
        annot_path,
        "GO_Bubble_down_CC_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top30_",
        projectname,
        "_200x250mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 250,
      width = 200,
      units = "mm",
      dpi = "retina"
    )
  } else
  {
    ego_down_CC_simple <- NULL
  }
  
  if(nrow(base::as.data.frame(ego_down_MF)) > 0)
  {
    ego_down_MF <- pairwise_termsim(ego_down_MF)
    ego_down_MF_simple <- clusterProfiler::simplify(ego_down_MF, 
                                                    cutoff = 0.7, 
                                                    by = "p.adjust", 
                                                    select_fun = min)
    # Bubble plot
    dotplot(ego_down_MF_simple, 
            x = "FoldEnrichment", 
            title = "Down - Molecular Function", 
            showCategory = 30)
    ggsave(
      paste0(
        annot_path,
        "GO_Bubble_down_MF_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top30_",
        projectname,
        "_200x250mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 250,
      width = 200,
      units = "mm",
      dpi = "retina"
    )
  } else
  {
    ego_down_MF_simple <- NULL
  }
  
  #### Write combined xlsx workbook ####
  # This combines the different GO tables both complete
  # and simplified
  
  wbook <- createWorkbook("GO_Enrich_Results")
  addWorksheet(wbook, "Up_BP")
  addWorksheet(wbook, "Up_CC")
  addWorksheet(wbook, "Up_MF")
  addWorksheet(wbook, "Down_BP")
  addWorksheet(wbook, "Down_CC")
  addWorksheet(wbook, "Down_MF")
  addWorksheet(wbook, "Up_BP_simple")
  addWorksheet(wbook, "Up_CC_simple")
  addWorksheet(wbook, "Up_MF_simple")
  addWorksheet(wbook, "Down_BP_simple")
  addWorksheet(wbook, "Down_CC_simple")
  addWorksheet(wbook, "Down_MF_simple")
  writeData(wbook, sheet = 1, ego_up_BP)
  writeData(wbook, sheet = 2, ego_up_CC)
  writeData(wbook, sheet = 3, ego_up_MF)
  writeData(wbook, sheet = 4, ego_down_BP)
  writeData(wbook, sheet = 5, ego_down_CC)
  writeData(wbook, sheet = 6, ego_down_MF)
  writeData(wbook, sheet = 7, ego_up_BP_simple)
  writeData(wbook, sheet = 8, ego_up_CC_simple)
  writeData(wbook, sheet = 9, ego_up_MF_simple)
  writeData(wbook, sheet = 10, ego_down_BP_simple)
  writeData(wbook, sheet = 11, ego_down_CC_simple)
  writeData(wbook, sheet = 12, ego_down_MF_simple)
  saveWorkbook(wbook,
               paste0(annot_path,"GO_overrep_", i, "_"
                      ,"pVal<0.05_qVal<0.2_Top30_", 
                      projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),
                      ".xlsx"),
               overwrite = TRUE)
  }

## GSEA analysis####
### Generate geneList 
# Should be list of SORTED fold changes named (names) by identifiers 
# (gene names etc.). I will base this on the candidate lists for individual
# comparisons.

dummy <- candidates_DS_K %>% filter(Gene.names != "")
geneList <- dummy$logFC #get numeric values
names(geneList) <- as.character(dummy$Gene.names) # get gene names
geneList <- sort(geneList, decreasing = TRUE) # sort list be numeric value

### Perform GSEA analysis
gseGO_BP <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              keyType = "SYMBOL",
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

heatplot(gseGO_BP, foldChange=geneList, showCategory=5) # because of duplicate 
#gene names coloring by fold change does not work, at least i think this is the cause

## KEGG-terms####
# First we generate the background lists (parameter: "universe") of all protein ids 
# followed by the individual set of differential enriched proteins for every 
# comparison (parameter: "gene")

### Over-representation analysis ####
#### Generate background list####
background_uniprot_ids <- psite_table$Leading.proteins %>%
  str_split(";", simplify = TRUE) %>%
  as.data.frame()  %>% 
  dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
background_uniprot_ids <- background_uniprot_ids$V1
background_uniprot_ids <- na.omit(background_uniprot_ids)

#### Generate DE-protein lists####
for (i in x) {
  limma_result <- topTable(fit3, coef = i, number = Inf)
  limma_result$ProtID <- row.names(limma_result)
  limma_result$MajorityProtID <- protein[rownames(limma_result), "Majority.protein.IDs"]
  
  up_uniprot_ids <- filter(limma_result,
                   logFC > 0.58 & adj.P.Val < 0.05 ) %>%
    dplyr::select(c("MajorityProtID"))
  up_uniprot_ids <- str_split(up_uniprot_ids$MajorityProtID, ";", simplify = TRUE) %>%
    as.data.frame() %>% 
    dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
  up_uniprot_ids <- up_uniprot_ids$V1
  up_uniprot_ids <- na.omit(up_uniprot_ids)
  
  down_uniprot_ids <- filter(limma_result,
                     logFC < -0.58 & adj.P.Val < 0.05 ) %>%
    dplyr::select(c("MajorityProtID"))
  down_uniprot_ids <- str_split(down_uniprot_ids$MajorityProtID, ";", simplify = TRUE) %>%
    as.data.frame() %>% 
    dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
  down_uniprot_ids <- down_uniprot_ids$V1
  down_uniprot_ids <- na.omit(down_uniprot_ids)
  
  #### Perform over-representation analyses####
  
  ego_down_KEGG <- enrichKEGG(
    gene          = down_uniprot_ids,
    universe      = background_uniprot_ids,
    keyType = "uniprot",
    organism       = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.2
  )
  
  ego_up_KEGG <- enrichKEGG(
    gene          = up_uniprot_ids,
    universe      = background_uniprot_ids,
    keyType = "uniprot",
    organism = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.2
  )
  
  #### Calculate fold enrichment and add it as an extra column####
  ego_down_KEGG@result <- ego_down_KEGG@result %>% 
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  
  ego_up_KEGG@result <- ego_up_KEGG@result %>% 
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  
  
  #### Plot data####
  
  # Bubble plot
  if(nrow(base::as.data.frame(ego_up_KEGG)) > 0)
  {  dotplot(ego_up_KEGG, 
             x = "FoldEnrichment", 
             title = "Up - KEGG", 
             showCategory = 30)
    ggsave(
      paste0(
        annot_path,
        "KEGG_Bubble_up_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top30_",
        projectname,
        "_200x250mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 250,
      width = 200,
      units = "mm",
      dpi = "retina"
    )
  } else
  {
    ego_up_KEGG <- NULL
  }
  
  # Bubble plot
  if(nrow(base::as.data.frame(ego_down_KEGG)) > 0)
  { dotplot(ego_down_KEGG, 
            x = "FoldEnrichment", 
            title = "Down - KEGG", 
            showCategory = 30)
    ggsave(
      paste0(
        annot_path,
        "KEGG_Bubble_down_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top30_",
        projectname,
        "_200x250mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 250,
      width = 200,
      units = "mm",
      dpi = "retina"
    )
  } else
  {
    ego_down_KEGG <- NULL
  }
  
  #### Write combined xlsx workbook ####
  # This combines the different GO tables both complete
  # and simplified
  
  wbook <- createWorkbook("GO_Enrich_Results")
  addWorksheet(wbook, "Up_KEGG")
  addWorksheet(wbook, "Down_KEGG")
  writeData(wbook, sheet = 1, ego_up_KEGG)
  writeData(wbook, sheet = 2, ego_down_KEGG)

  saveWorkbook(wbook,
               paste0(annot_path,"KEGG_overrep_", i, "_"
                      ,"pVal<0.05_qVal<0.2_Top30_", 
                      projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),
                      ".xlsx"),
               overwrite = TRUE)
}
