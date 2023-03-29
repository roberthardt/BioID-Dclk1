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
    "ggpubr",
    "VennDiagram"
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


# Set output paths for plots and tables
QC_path <- "P-Sites/Output/QC_plots/"

# Select file
psites <-
  read.table(
    file = "P-Sites/Data/txt_03_BioID-DCLK1_200408_Phospho (STY)Sites.txt",
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
dat <- dplyr::select(p_filter, contains(c(
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


# Valid values filtering#######################################################
# Count number of NaN per condition.
dat_log$na_count_DS_L = apply(dat_log, 1, function(x)
  sum(is.na(x[c(1:3)])))
dat_log$na_count_DSK_L = apply(dat_log, 1, function(x)
  sum(is.na(x[c(4:6)])))
dat_log$na_count_K_L = apply(dat_log, 1, function(x)
  sum(is.na(x[c(7:9)])))
dat_log$na_count_SK_L = apply(dat_log, 1, function(x)
  sum(is.na(x[c(10:12)])))


# Filter table for valid values. Here we want a minimum of 2 values 
# in any Light sample group.
dat_log_filter <- dat_log %>%
  dplyr::filter(if_any(contains("na_count"),
                ~ . <= 1)) %>%
  dplyr::select(c(1:24))

# Append data to sites count table
count_tbl[nrow(count_tbl)+1,] <- c("2_valid_L_val_anygroup",
                                   nrow(dat_log_filter),
                                   "ClassI sites quantified in light channel for at least 2 times in any sample group")

################################################################################
# Imputation####################################################################
################################################################################

# Impute missing signals by the lowest value in all channels (min)
# First store position info of imputed values in a separate data frame
na_pos <- is.na(dat_log_filter) %>%
  as.data.frame()

# Set seed
set.seed(1234)

# Perform imputation
dat_log_impute <- dat_log_filter %>%
  as.matrix() %>%
  MsCoreUtils::impute_matrix(method = "min") %>%
  as.data.frame()

################################################################################
# QC Plots######################################################################
################################################################################

# Plot density plots for data after imputation
# Generate some basic qc plots (box, density) in ggplot. 
## Reformatting####
### Therefore convert to long format first
dat_log_impute_long <- dat_log_impute %>%
  rownames_to_column(var = "SiteID") %>%
  pivot_longer(
    cols = -SiteID,
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


#Pivot also the NA-Positions to the long format
na_pos_long <- na_pos %>%
  rownames_to_column(var = "SiteID") %>%
  pivot_longer(
    cols = -SiteID,
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

#Merge expression data frame and na_position frame
dat_log_impute_merge <- merge(x = dat_log_impute_long,
                              y = na_pos_long[, c("sampleID", "SiteID", "imputed")],
                              by = c("sampleID","SiteID"))

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
  facet_grid(label ~ sampleID, scales = "free") +
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
    "P-sites_Density_log2_intensities+impute_L+H_400x150_",
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
    "P-sites_Density_log2_intensities+impute_L+H_400x150_",
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

################################################################################
# Ratio generation and further filtering########################################
################################################################################

## Generate L/H ratios####
df_ratios <- dat_log_impute[, 1:12] - dat_log_impute[, 13:24]

## Rename colnames####
colnames(df_ratios) <- c("log2_ratio_LH_DS_01",
                         "log2_ratio_LH_DS_02",
                         "log2_ratio_LH_DS_03",
                         "log2_ratio_LH_DSK_01",
                         "log2_ratio_LH_DSK_02",
                         "log2_ratio_LH_DSK_03",
                         "log2_ratio_LH_K_01",
                         "log2_ratio_LH_K_02",
                         "log2_ratio_LH_K_03",
                         "log2_ratio_LH_SK_01",
                         "log2_ratio_LH_SK_02",
                         "log2_ratio_LH_SK_03")

## Filter for enriched proteins ####
### Count number of log2 L/H ratios >1 per condition.####
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

### Generate separate tables of enriched proteins per construct ####
# Also add the corresponding log2 abundances before imputation

#### DSK ####
DSK_enr <- df_ratios %>%
  filter(upreg_count_DSK >=2) %>%
  dplyr::select(contains("log2_ratio_LH_DSK_")) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_log_filter %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("DSK_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>% 
  left_join(x = ., 
            y = dat_log %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("NA_count_DSK_L"))) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>% 
  left_join(x = ., 
            y = psites %>%
              dplyr::select(c("Site", "Leading.proteins","Gene.names", "Amino.acid", 
                              "Position", "Fasta.headers")) %>%
              dplyr::rename(SiteID = Site),
            by = "SiteID")

#### DS####
DS_enr <- df_ratios %>%
  filter(upreg_count_DS >=2) %>%
  dplyr::select(contains("log2_ratio_LH_DS_")) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_log_filter %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("DS_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>% 
  left_join(x = ., 
            y = dat_log %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("NA_count_DS_L"))) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>% 
  left_join(x = ., 
            y = psites %>%
              dplyr::select(c("Site", "Leading.proteins","Gene.names", "Amino.acid", 
                              "Position", "Fasta.headers")) %>%
              dplyr::rename(SiteID = Site),
            by = "SiteID")

#### SK ####
SK_enr <- df_ratios %>%
  filter(upreg_count_SK >=2) %>%
  dplyr::select(contains("log2_ratio_LH_SK_")) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_log_filter %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("SK_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>% 
  left_join(x = ., 
            y = dat_log %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("NA_count_SK_L"))) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID")  %>% 
  left_join(x = ., 
            y = psites %>%
              dplyr::select(c("Site", "Leading.proteins","Gene.names", "Amino.acid", 
                              "Position", "Fasta.headers")) %>%
              dplyr::rename(SiteID = Site),
            by = "SiteID")

#### K ####
K_enr <- df_ratios %>%
  filter(upreg_count_K >=2) %>%
  dplyr::select(contains("log2_ratio_LH_K_")) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_log_filter %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("K_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>% 
  left_join(x = ., 
            y = dat_log %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("NA_count_K_L"))) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>% 
  left_join(x = ., 
            y = psites %>%
              dplyr::select(c("Site", "Leading.proteins","Gene.names", "Amino.acid", 
                              "Position", "Fasta.headers")) %>%
              dplyr::rename(SiteID = Site),
            by = "SiteID")


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
                                   "ClassI sites having at least 2 log2 L/H ratios >1 in any sample group")

################################################################################
# Write biotin-enriched output table for paper #################################
################################################################################

## Create output subfolder
if (dir.exists("P-Sites/Output/Biotin_enriched")) {
  
  cat("The folder already exists")
  
} else {
  
  dir.create("P-Sites/Output/Biotin_enriched",
             recursive = TRUE)
  
}

## Abundance table####
### Identified ClassI sites####
# Generate p-site table with abundance values from log2 abundance matrix
# Add log2 abundances
psite_id_table <-
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
    y = dat_log %>%
      dplyr::select(contains(c("01", "02", "03"))) %>%
      as.data.frame() %>%
      rename_with(.cols = everything(), 
                  function(x){paste0("log2_intensity_", x)}) %>%
      rownames_to_column(var = "Site"),
    by = "Site"
  )


### Quantified ClassI sites####
# Generate p-site table with abundance values from log2_filter abundance matrix
# Add log2_filter abundances. These are sites quantified with 2 valid values 
# in at least 1 construct, We also add the imputed abundances and the resulting 
# L/H-ratios
psite_quant_table <-
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
    y = dat_log_filter %>%
      as.data.frame() %>%
      rename_with(.cols = everything(), 
                  function(x){paste0("log2_intensity_", x)}) %>%
      rownames_to_column(var = "Site"),
    by = "Site"
  ) %>%
  left_join(x = ., 
            y = dat_log_impute %>%
              as.data.frame() %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_imputed_intensity_", x)}) %>%
              rownames_to_column(var = "Site"),
            by = "Site") %>%
  left_join(x = .,
            y = na_pos %>% 
              rename_with(.cols = everything(), 
                          function(x){paste0("is_imputed_", x)}) %>%
              rownames_to_column(var = "Site"),
            by = "Site"
  ) %>%
  left_join(x = ., 
            y = df_ratios %>%
              as.data.frame() %>%
              dplyr::select(contains("ratio")) %>%
              rownames_to_column(var = "Site"),
            by = "Site")

### Enriched ClassI sites####
#Generate p-site table with abundance values from log2_enriched abundance matrix 
# also add the L/H-ratios
psite_enr_table <-
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
    y = df_ratios_filter %>%
      as.data.frame() %>%
      dplyr::select(contains("ratio")) %>%
      rownames_to_column(var = "Site"),
    by = "Site") %>%
  left_join(x = ., 
            y = dat_log_enriched %>%
              as.data.frame() %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_intensity_", x)}) %>%
              rownames_to_column(var = "Site"),
            by = "Site")


## Combined biotin-enriched analysis output workbook####
# Write xlsx workbook which combines psite_tables and the candidate tables
# of biotin-enriched psites
wbook <- createWorkbook("P-sites_Biotin-Enriched_Results")
addWorksheet(wbook, "Identified_P-sites")
addWorksheet(wbook, "Quantified_P-sites")
addWorksheet(wbook, "Enriched_P-sites")
addWorksheet(wbook, "enr_DSK")
addWorksheet(wbook, "enr_DS")
addWorksheet(wbook, "enr_SK")
addWorksheet(wbook, "enr_K")
writeData(wbook, sheet = 1, psite_id_table)
writeData(wbook, sheet = 2, psite_quant_table)
writeData(wbook, sheet = 3, psite_enr_table)
writeData(wbook, sheet = 4, DSK_enr)
writeData(wbook, sheet = 5, DS_enr)
writeData(wbook, sheet = 6, SK_enr)
writeData(wbook, sheet = 7, K_enr)
saveWorkbook(wbook,
             paste0("P-Sites/Output/Biotin_enriched/", "Biotin-enriched_results_",
                    projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".xlsx"),
             overwrite = TRUE)

#############################################################################
# Annotation enrichment analysis of Biotin-enriched proteins ##################
#############################################################################

## Load packages ####
library(clusterProfiler) #load clusterProfiler package
library(org.Mm.eg.db) # load annotation database for Mus musculus
library(DOSE) # load DOSE package for parsing results
library(enrichplot) # needed for GO plots

## Create output subfolder
if (dir.exists("P-Sites/Output/Biotin_enriched/Annotation_enrichment")) {
  
  cat("The folder already exists")
  
} else {
  
  dir.create("P-Sites/Output/Biotin_enriched/Annotation_enrichment",
             recursive = TRUE)
  
}

## Set relative output path####
Biotin_annot_path <- "P-Sites/Output/Biotin_enriched/Annotation_enrichment/"

## GO-terms####
# Since we want to see what kind of proteins are interacting with the different
# Dclk1-constructs, we compare the enriched proteins against the global background

### Over-representation analysis ####
#### Generate DE-protein lists####
y <- c("DSK_enr", "DS_enr", "SK_enr", "K_enr")

for (i in y) {
  my_ids <- dplyr::select(get(i), c("Gene.names"))
  my_ids <- str_split(my_ids$Gene.names, ";", simplify = TRUE) %>%
    as.data.frame() %>% 
    dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
  my_ids <- my_ids$V1
  my_ids <- na.omit(my_ids)
  my_ids <- unique(my_ids)

  
  #### Perform over-representation analyses####
  ego_BP <- enrichGO(
    gene          = my_ids,
    keyType = "SYMBOL",
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.2,
    readable      = FALSE
  )
  
  ego_CC <- enrichGO(
    gene          = my_ids,
    keyType = "SYMBOL",
    OrgDb         = org.Mm.eg.db,
    ont           = "CC",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.2,
    readable      = FALSE
  )
  ego_MF <- enrichGO(
    gene          = my_ids,
    keyType = "SYMBOL",
    OrgDb         = org.Mm.eg.db,
    ont           = "MF",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.2,
    readable      = FALSE
  )
  
  
  #### Calculate fold enrichment and add it as an extra column####
  ego_BP@result <- ego_BP@result %>% 
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  
  ego_CC@result <- ego_CC@result %>% 
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  
  ego_MF@result <- ego_MF@result %>% 
    mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  

  #### Simplify terms and plot####
  # GO-Terms are simplified by removing redundant terms (semantic similarity)
  # Afterwards plots are generated if enriched terms are found
  
  if(nrow(base::as.data.frame(ego_BP@result)) > 0)
  {
    ego_BP <- enrichplot::pairwise_termsim(ego_BP)
    ego_BP_simple <- clusterProfiler::simplify(ego_BP, 
                                                  cutoff = 0.7, 
                                                  by = "p.adjust", 
                                                  select_fun = min)
    # Bubble plot
    dotplot(ego_BP_simple, 
            x = "FoldEnrichment", 
            showCategory = 25,
            title = paste("Biological Process", i, sep = " - "))+
      theme(
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold")
      )
    ggsave(
      paste0(
        Biotin_annot_path,
        "GO_Bubble_BP_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top25_",
        projectname,
        "_200x260mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 260,
      width = 200,
      units = "mm",
      dpi = "retina"
    )
  } else
  {
    ego_BP_simple <- NULL
  }
  
  if(nrow(base::as.data.frame(ego_CC)) > 0)
  {
    ego_CC <- pairwise_termsim(ego_CC)
    ego_CC_simple <- clusterProfiler::simplify(ego_CC, 
                                                  cutoff = 0.7, 
                                                  by = "p.adjust", 
                                                  select_fun = min)
    # Bubble plot
    dotplot(ego_CC_simple, 
            x = "FoldEnrichment", 
            title = paste("Cellular Component", i, sep = " - "), 
            showCategory = 25)+
      theme(
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))
    
    ggsave(
      paste0(
        Biotin_annot_path,
        "GO_Bubble_CC_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top25_",
        projectname,
        "_200x260mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 260,
      width = 200,
      units = "mm",
      dpi = "retina"
    )
  } else
  {
    ego_CC_simple <- NULL
  }
  
  if(nrow(base::as.data.frame(ego_MF)) > 0)
  {
    ego_MF <- pairwise_termsim(ego_MF)
    ego_MF_simple <- clusterProfiler::simplify(ego_MF, 
                                                  cutoff = 0.7, 
                                                  by = "p.adjust", 
                                                  select_fun = min)
    # Bubble plot
    dotplot(ego_MF_simple, 
            x = "FoldEnrichment", 
            title = paste("Molecular Function", i, sep = " - "), 
            showCategory = 25)+
      theme(
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))
    
    ggsave(
      paste0(
        Biotin_annot_path,
        "GO_Bubble_MF_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top25_",
        projectname,
        "_200x260mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 260,
      width = 200,
      units = "mm",
      dpi = "retina"
    )
  } else
  {
    ego_MF_simple <- NULL
  }
  
  
  
  #### Write combined xlsx workbook ####
  # This combines the different GO tables both complete
  # and simplified
  
  wbook <- createWorkbook("GO_Enrich_Results")
  addWorksheet(wbook, "BP")
  addWorksheet(wbook, "CC")
  addWorksheet(wbook, "MF")
  addWorksheet(wbook, "BP_simple")
  addWorksheet(wbook, "CC_simple")
  addWorksheet(wbook, "MF_simple")
  writeData(wbook, sheet = 1, ego_BP)
  writeData(wbook, sheet = 2, ego_CC)
  writeData(wbook, sheet = 3, ego_MF)
  writeData(wbook, sheet = 4, ego_BP_simple)
  writeData(wbook, sheet = 5, ego_CC_simple)
  writeData(wbook, sheet = 6, ego_MF_simple)
  saveWorkbook(wbook,
               paste0(Biotin_annot_path,"GO_overrep_", i, "_"
                      ,"pVal<0.05_qVal<0.2_", 
                      projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),
                      ".xlsx"),
               overwrite = TRUE)
}

################################################################################
# Venn of Biotin-enriched vs. Dimethyl dataset##################################
################################################################################
# We want to compare the Biotin-enriched phosphosites against the ones regulated
# in the dimethyl dataset day7/0 and day 14/0.

## Extract gene names for each regulation ####
### biotin-enriched = this dataset
biotin_enr_ids <- psite_enr_table %>%
  dplyr::select(c("Gene.names", "Amino.acid", "Position")) %>% 
  unite(col = Site, 
        c(Gene.names, Amino.acid), 
        remove = T) %>% unite(col = Site,
                              c(Site,Position),
                              sep = "",
                              remove = T) %>%
  separate_wider_delim(cols = everything(), 
                       delim = ";",
                       names = c("Site"),
                       too_many = "drop") %>% 
  filter(Site != "")

### day7 vs 0 - dimethyl dataset
DML_7_0_all <- read.xlsx(xlsxFile = 
                           "P-Sites/Data/Table_S3_Neurospheres_Phosphoproteome_processed_data.xlsx", 
                         sheet = 4) 
DML_7_0_reg <- DML_7_0_all %>% 
  dplyr::filter (Regulated == "+") %>%
  dplyr::select(c("Gene.names", "Amino.acid", "Position")) %>%
  unite(col = Site, 
        c(Gene.names, Amino.acid), 
        remove = T) %>% unite(col = Site,
                                  c(Site,Position),
                                  sep = "",
                                  remove = T) %>%
  distinct() %>%
  separate_wider_delim(cols = everything(), 
                       delim = ";",
                       names = c("Site"),
                       too_many = "drop") %>% 
  filter(Site != "")

### day14 vs 0 - dimethyl dataset
DML_14_0_all <- read.xlsx(xlsxFile = 
                            "P-Sites/Data/Table_S3_Neurospheres_Phosphoproteome_processed_data.xlsx", 
                          sheet = 5) 
DML_14_0_reg <- DML_14_0_all %>% 
  dplyr::filter (Regulated == "+") %>%
  dplyr::select(c("Gene.names", "Amino.acid", "Position")) %>%
  unite(col = Site, 
        c(Gene.names, Amino.acid), 
        remove = T) %>% unite(col = Site,
                              c(Site,Position),
                              sep = "",
                              remove = T) %>%
  distinct() %>%
  separate_wider_delim(cols = everything(), 
                       delim = ";",
                       names = c("Site"),
                       too_many = "drop") %>% 
  filter(Site != "")

## Create lists for each Venn to be ploted ####
venn_a <- list( all_BioID = biotin_enr_ids$Site,
                 day7_0 = DML_7_0_reg$Site,
                 day14_0 = DML_14_0_reg$Site)

## Calculate Venn overlaps & write to file
venn_a_counts <- get.venn.partitions(venn_a)
write.xlsx(x = venn_a_counts, 
           file = paste0("P-Sites/Output/Venn_diagrams/Venn-table_BioID-Biotin-riched-vs-DML_",projectname,
                         "_",format(Sys.time(), '%Y%m%d_%H%M%S'),".xlsx"))

## Plot Venn diagrams ####
venn.diagram(venn_a,
             filename = paste0("Venn_BioID-Biotin-riched-vs-DML_",
                               projectname, "_",format(Sys.time(), 
                                                       '%Y%m%d_%H%M%S'),".png"),
             category.names = c(all_BioID = "All BioID",
                                day7_0 = "Day 7 vs. day 0",
                                day14_0 = "Day 14 vs. day 0"),
             print.mode = c("raw", "percent"),
             sigdigs = 1,
             fill = c("#EFC000FF", "#5CB85CFF", "#9632B8FF"),
             fontfamily = "sans",
             cat.fontfamily = "sans",
             cat.fontface = "bold",
             cat.col = c("#EFC000FF", "#5CB85CFF", "#9632B8FF"),
             cat.cex = 1.4,
             cat.just = list(a = c(0.5, -1),
                             b = c(0.5, -1),
                             c = c(0.5, 1)),
             imagetype = "png",
             margin = 0.1)


# QC Plots ---------------------------------------------------------------------
# Plot light log2 intensities in boxplot and density plots of the 
# enriched proteins
## Reformatting#####
### Convert to long format
dat_log_enriched_long <- dat_log_enriched %>%
  rownames_to_column(var = "siteID") %>%
  pivot_longer(
    cols = -siteID,
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

################################################################################
# Normalization ################################################################
################################################################################

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
  write_tsv(file = "P-Sites/Output/normalizerDE_eval/designNormalyzerDE2.tsv")

# Generate input data & write to file.
# Input data has to be "de-logarithmized" before
input_data <- 2 ^ dat_log_enriched  %>%
  mutate(ProtID = row.names(dat_log_enriched),
         .before = 1) %>%
  write_tsv(file = "P-Sites/Output/normalizerDE_eval/inputNormalyzerDE.tsv")

# Run normalyzerDE
NormalyzerDE::normalyzer(
  jobName = "normalizerDE_eval",
  designPath = "P-Sites/Output/normalizerDE_eval/designNormalyzerDE.tsv",
  dataPath = "P-Sites/Output/normalizerDE_eval/inputNormalyzerDE.tsv",
  outputDir = "P-Sites/Output/normalizerDE_eval/"
)

# Result: ==> basically every method is better then simple log2 or mean. But RLR,
# cyclicLoess and vsn perform best. Cyclic loess leads to lower Intragroup
# Pearson & Spearman correlations. I opt for RLR.

## Perform normalization ####
dat_norm <-
  performGlobalRLRNormalization(as.matrix(dat_log_enriched), 
                                noLogTransform = TRUE) %>%
  as.data.frame()

################################################################################
# QC Plots #####################################################################
################################################################################

# Plot normalized data in Boxplot, Density plot and PCA

## Reformatting ####
# Convert to long format
dat_norm_long <- dat_norm %>%
  as.data.frame() %>%
  mutate(siteID = row.names(dat_norm)) %>%
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

###############################################################################
# Imputation-################################################################
# We impute missing values of light intensities to get more robust statistics
###############################################################################

# Impute missing signals by normal distribution
#First store position info of imputed values in a separate data frame
na_pos_2 <- is.na(dat_norm) %>%
  as.data.frame()

# Set seed to always get the same imputations when runnig the script 
set.seed(1234)

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
  rownames_to_column(var = "SiteID") %>%
  pivot_longer(
    cols = -SiteID,
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


#Pivot also the NA-Positions to the long format
na_pos_long_2 <- na_pos_2 %>%
  rownames_to_column(var = "SiteID") %>%
  pivot_longer(
    cols = -SiteID,
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


#Merge expression data frame and na_position frame
dat_norm_impute_merge <- merge(x = dat_norm_impute_long,
                              y = na_pos_long_2[, c("sampleID", "SiteID", "imputed")],
                              by = c("sampleID", "SiteID"))

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

#### all constructs ####

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

#### w/o DS ####

#1. Prepare data for PCA

#Remove DS from table
mydata_noDS <- dat_norm_impute %>%
  dplyr::select(contains(c("DSK", "SK", "K")))

#Remove NaN
mydata_noNaN <- na.omit(mydata_noDS)

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
  scale_color_manual(values = c("#1B9E77","#7570B3","#E7298A"),
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
         "P-sites_PCA_wo-DS_",
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
         "P-sites_PCA_wo-DS_",
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

#### all constructs ####
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
  annotation_colors = heatcolors,
  fontsize = 16,
  angle_col = 45 
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

#### w/o DS ####
#Generate annotation table
heatanno <- data.frame(sample = colnames(dat_norm_impute))
heatanno$Construct <- str_extract(heatanno$sample, "[^_]+")
rownames(heatanno) <- colnames(dat_norm_impute)
heatanno <- subset (heatanno, select = -sample)
heatanno$Construct <- factor(heatanno$Construct, levels = c("DSK", "DS", "SK", "K"))


# List with colors for each annotation.
heatcolors <- list(Construct = c("#1B9E77","#7570B3","#E7298A"))
names(heatcolors$Construct) <- c("DSK", "SK", "K")

plotheat <- pheatmap(
  dat_norm_impute %>% dplyr::select((contains(c("DSK", "SK", "K")))),
  clustering_distance_cols = "manhattan",
  show_rownames = FALSE,
  cutree_cols = 3,
  color = colorRampPalette(c("navy", "white", "red"))(50),
  na_col = "grey39",
  annotation_col = heatanno,
  annotation_colors = heatcolors,
  fontsize = 16,
  angle_col = 45 
)

png(
  filename = paste0(QC_path,
                    "P-sites_Heatmap_wo-DS_",
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

DE_path <- "P-Sites/Output/DE_analysis/" # Set path for output data

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
    ylim = c(0, max(-log10(limma_result$P.Value) + 1)),
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
      "_200x200mm_",
      format(Sys.time(), '%Y%m%d_%H%M%S'),
      ".png"
    ),
    height = 200,
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
      size = 5,
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
      "P-sites_MDplot_",
      gsub(" ", "", i),
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

################################################################################
# Write DE-output table for paper ##############################################
################################################################################

## Candidate tables####
all_DSK_DS <- limma::topTable(fit3, coef = "DSK-DS", number = Inf)
candidates_DSK_DS <- all_DSK_DS %>%
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(`logFC` > 0.58 & `adj.P.Val` < 0.05 |
                  `logFC` < -0.58 & `adj.P.Val` < 0.05 ) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_norm_impute %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("DSK_","DS_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_norm_imputed_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>%
  left_join(x = ., 
            y = psites %>%
              dplyr::select(c("Site", "Gene.names", "Amino.acid", "Position", "Fasta.headers")) %>%
              dplyr::rename(SiteID = Site))

all_DSK_SK <- topTable(fit3, coef = "DSK-SK", number = Inf)
candidates_DSK_SK <- all_DSK_SK %>%
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(`logFC` > 0.58 & `adj.P.Val` < 0.05 |
                  `logFC` < -0.58 & `adj.P.Val` < 0.05 ) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_norm_impute %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("DSK_","SK_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_norm_imputed_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>%
  left_join(x = ., 
            y = psites %>%
              dplyr::select(c("Site", "Gene.names", "Amino.acid", "Position", "Fasta.headers")) %>%
              dplyr::rename(SiteID = Site))

all_DSK_K <- topTable(fit3, coef = "DSK-K", number = Inf) 
candidates_DSK_K <- all_DSK_K %>%
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(`logFC` > 0.58 & `adj.P.Val` < 0.05 |
                  `logFC` < -0.58 & `adj.P.Val` < 0.05 ) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_norm_impute %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("DSK_","K_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_norm_imputed_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>%
  left_join(x = ., 
            y = psites %>%
              dplyr::select(c("Site", "Gene.names", "Amino.acid", "Position", "Fasta.headers")) %>%
              dplyr::rename(SiteID = Site))

all_DS_SK <- topTable(fit3, coef = "DS-SK", number = Inf)
candidates_DS_SK <- all_DS_SK %>%
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(`logFC` > 0.58 & `adj.P.Val` < 0.05 |
                  `logFC` < -0.58 & `adj.P.Val` < 0.05 ) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_norm_impute %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("DS_","SK_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_norm_imputed_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>%
  left_join(x = ., 
            y = psites %>%
              dplyr::select(c("Site", "Gene.names", "Amino.acid", "Position", "Fasta.headers")) %>%
              dplyr::rename(SiteID = Site))

all_DS_K <- topTable(fit3, coef = "DS-K", number = Inf)
candidates_DS_K <- all_DS_K %>%
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(`logFC` > 0.58 & `adj.P.Val` < 0.05 |
                  `logFC` < -0.58 & `adj.P.Val` < 0.05 ) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_norm_impute %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("DS_","K_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_norm_imputed_intensity_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>%
  left_join(x = ., 
            y = psites %>%
              dplyr::select(c("Site", "Gene.names", "Amino.acid", "Position", "Fasta.headers")) %>%
              dplyr::rename(SiteID = Site))

all_SK_K <- topTable(fit3, coef = "SK-K", number = Inf)
candidates_SK_K <- all_SK_K %>%
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
  dplyr::filter(`logFC` > 0.58 & `adj.P.Val` < 0.05 |
                  `logFC` < -0.58 & `adj.P.Val` < 0.05 ) %>%
  rownames_to_column(var ="SiteID") %>% 
  left_join(x = ., 
            y = dat_norm_impute %>%
              as.data.frame() %>%
              dplyr::select(starts_with(c("SK_","K_"))) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0("log2_norm_intensity_imputed_", x)}) %>%
              rownames_to_column(var = "SiteID"),
            by = "SiteID") %>%
  left_join(x = ., 
            y = psites %>%
              dplyr::select(c("Site", "Gene.names", "Amino.acid", "Position", "Fasta.headers")) %>%
              dplyr::rename(SiteID = Site))

## Abundance table####
#Generate protein table with logFCs, p.values and abundance values from 
# normalized matrix and norm+imputed matrix
# First add the DE results
psite_norm_table <-
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
                      "PEP")) %>%
      dplyr::rename(SiteID = Site),
    y = all_DSK_DS %>%
      as.data.frame() %>%
      dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
      rename_with(.cols = everything(), 
                  function(x){paste0(x, "_DSK_DS")}) %>%
      rownames_to_column(var = "SiteID") %>%
      dplyr::mutate(significant_DSK_DS = SiteID %in% 
                      candidates_DSK_DS$SiteID),
    by = "SiteID") %>%
  left_join(x = .,
            y = all_DSK_SK %>%
              as.data.frame() %>%
              dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0(x, "_DSK_SK")}) %>%
              rownames_to_column(var = "SiteID") %>%
              dplyr::mutate(significant_DSK_SK = SiteID %in% 
                              candidates_DSK_SK$SiteID),
            by = "SiteID") %>%
  left_join(x = .,
            y = all_DSK_K %>%
              as.data.frame() %>%
              dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0(x, "_DSK_K")}) %>%
              rownames_to_column(var = "SiteID") %>%
              dplyr::mutate(significant_DSK_K = SiteID %in% 
                              candidates_DSK_K$SiteID),
            by = "SiteID") %>%
  left_join(x = .,
            y = all_DS_SK %>%
              as.data.frame() %>%
              dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0(x, "_DS_SK")}) %>%
              rownames_to_column(var = "SiteID") %>%
              dplyr::mutate(significant_DS_SK = SiteID %in% 
                              candidates_DS_SK$SiteID),
            by = "SiteID") %>%
  left_join(x = .,
            y = all_DS_K %>%
              as.data.frame() %>%
              dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0(x, "_DS_K")}) %>%
              rownames_to_column(var = "SiteID") %>%
              dplyr::mutate(significant_DS_K = SiteID %in% 
                              candidates_DS_K$SiteID),
            by = "SiteID") %>%
  left_join(x = .,
            y = all_SK_K %>%
              as.data.frame() %>%
              dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
              rename_with(.cols = everything(), 
                          function(x){paste0(x, "_SK_K")}) %>%
              rownames_to_column(var = "SiteID") %>%
              dplyr::mutate(significant_SK_K = SiteID %in% 
                              candidates_SK_K$SiteID),
            by = "SiteID")

# Second add the normalized abundances
psite_norm_table <-
  left_join(
    x = psite_norm_table,
    y = dat_norm %>%
      as.data.frame() %>%
      rename_with(.cols = everything(), 
                  function(x){paste0("log2_norm_intensity_", x)}) %>%
      rownames_to_column(var = "SiteID"),
    by = "SiteID"
  )

#Then also add the imputed abundances
psite_norm_table <-
  left_join(
    x = psite_norm_table,
    y = dat_norm_impute %>%
      as.data.frame() %>%
      rename_with(.cols = everything(), 
                  function(x){paste0("log2_norm_imputed_intensity_", x)}) %>%
      rownames_to_column(var = "SiteID"),
    by = "SiteID"
  )

#Add indicator columns marking imputed values
psite_norm_table <-
  left_join(
    x = psite_norm_table,
    y = na_pos_2 %>% 
      rename_with(.cols = everything(), 
                  function(x){paste0("is_imputed_", x)}) %>%
      rownames_to_column(var = "SiteID"),
    by = "SiteID"
  )


## Combined DE-analysis output workbook####
# Write xlsx workbook which combines psite_norm_table and the candidate tables
wbook <- createWorkbook("ClassI_P-site_Results")
addWorksheet(wbook, "Enriched_P-Sites_norm")
addWorksheet(wbook, "Candidates_DSKvsDS")
addWorksheet(wbook, "Candidates_DSKvsSK")
addWorksheet(wbook, "Candidates_DSKvsK")
addWorksheet(wbook, "Candidates_DSvsSK")
addWorksheet(wbook, "Candidates_DSvsK")
addWorksheet(wbook, "Candidates_SKvsK")
writeData(wbook, sheet = 1, psite_norm_table)
writeData(wbook, sheet = 2, candidates_DSK_DS)
writeData(wbook, sheet = 3, candidates_DSK_SK)
writeData(wbook, sheet = 4, candidates_DSK_K)
writeData(wbook, sheet = 5, candidates_DS_SK)
writeData(wbook, sheet = 6, candidates_DS_K)
writeData(wbook, sheet = 7, candidates_SK_K)
saveWorkbook(wbook,
             paste0(DE_path, "Differential-expression_results_",
                    projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".xlsx"),
             overwrite = TRUE)


# Dclk1 detailed psite analysis ####----------------------------------------

## Set relative output path####
Dclk1_psites_path <- "P-Sites/Output/Dclk1_psite_analysis/"

## Filter for Dclk1 psites #####
Dclk1_sites <- psite_norm_table %>% 
  filter(grepl("Q9JLM8",Leading.proteins))

## Create new site label ####
Dclk1_sites <- Dclk1_sites %>% 
  unite(col = Site2, c(Amino.acid, Position),
        remove = FALSE)

## Select only site annotation and norm_intensity columns ####
Dclk1_sites_ma <- dplyr::select(Dclk1_sites, 
                         c(Position, 
                           Site2, 
                           contains("norm_intensity"))) %>% 
  dplyr:::arrange(Position)


## Remove unnecessary columns ####
Dclk1_sites_ma <- Dclk1_sites_ma %>% 
  dplyr::select(contains(c("norm_intensity","Site2")))

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
                            "SK_L_03",
                            "Site2")

### Arrange columns in specific order ####
Dclk1_sites_ma <- Dclk1_sites_ma %>%
  dplyr::select("DSK_L_01",
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
         "K_L_03",
         "Site2")

### Remove DS samples####
Dclk1_sites_ma_2 <- Dclk1_sites_ma %>%
  dplyr::select("DSK_L_01",
         "DSK_L_02",
         "DSK_L_03",
         "SK_L_01",
         "SK_L_02",
         "SK_L_03",
         "K_L_01",
         "K_L_02",
         "K_L_03",
         "Site2") 

### Remove rows with less then 2 values####
Dclk1_sites_ma_2 <- Dclk1_sites_ma_2[rowSums(is.na(Dclk1_sites_ma_2[1:9])) < 8, ]

### Generate Heatmap####
#### Generate annotation table####
heatanno <- data.frame(sample = colnames(Dclk1_sites_ma[1:12]))
heatanno$Construct <- rep(c("DSK", "DS", "SK", "K"),
                          each = 3)
rownames(heatanno) <- colnames(Dclk1_sites_ma[1:12])
heatanno <- subset (heatanno, select = -sample)
heatanno$Construct <- factor(heatanno$Construct, levels = c("DSK", "DS", "SK", "K"))


#### List with colors for each annotation####
heatcolors <- list(Construct = brewer.pal(4, "Dark2"))
names(heatcolors$Construct) <- c("DSK", "DS", "SK", "K")

#### List with colors without group DS####
heatcolors2 <- list(Construct = heatcolors$Construct[-2])


## Set rownames needed for heatmap ####
Dclk1_sites_ma %<>% as.data.frame()
row.names(Dclk1_sites_ma) <- Dclk1_sites_ma$Site2
Dclk1_sites_ma_2 %<>% as.data.frame()
row.names(Dclk1_sites_ma_2) <- Dclk1_sites_ma_2$Site2

#### Plot & Save heatmap - all#####
plotheat <- pheatmap(
  Dclk1_sites_ma[1:12],
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
  filename = paste0(Dclk1_psites_path,
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
  filename = paste0(Dclk1_psites_path,
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
  t(Dclk1_sites_ma_2[1:9]),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cutree_cols = 4,
  color = colorRampPalette(c("navy", "white", "red"))(50),
  na_col = "grey39",
  scale = "column",
  annotation_row = heatanno,
  annotation_colors = heatcolors2,
  fontsize = 16,
  angle_col = 45,
  cellwidth = 30,
  cellheight = 30
)

png(
  filename = paste0(Dclk1_psites_path,
                    "Dclk1_Sitemap_",
                    "enriched_norm_log2_intensities_L_z-scaled_woDS",
                    "_12x5in_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".png"
  ),
  width = 12 ,
  height = 5,
  units = "in",
  res = 300
)
plotheat
dev.off()

svg(
  filename = paste0(Dclk1_psites_path,
                    "Dclk1_Sitemap_",
                    "enriched_norm_log2_intensities_z-scaled_woDS",
                    "_8x6in_",
                    projectname,
                    "_",
                    format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".svg"
  ),
  width = 12 ,
  height = 5
)
plotheat
dev.off()


# Regulated psites heatmaps####-----------------------------------------------
## Filter prot_table for all regulated psites w/o DS####
psite_regulated <- psite_norm_table %>% 
  filter(SiteID %in% 
           c(candidates_DSK_K$SiteID,
             candidates_DSK_SK$SiteID,
             candidates_SK_K$SiteID))

# Append data to sites count table
count_tbl[nrow(count_tbl)+1,] <- c("regulated",
                                   nrow(psite_regulated),
                                   "ClassI sites showing significant regulation in any pairwise group comparison")

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

### Remove DS data ####
psite_reg_long %<>% filter(Sample_type != "DS")

### Generate summary table ####
psite_reg_count <- psite_reg_long %>%
  dplyr::select(Sample_type, Rep, log2_norm_intensity) %>%
  dplyr::group_by(Sample_type, Rep) %>%
  dplyr::summarise(Count = sum(!is.na(log2_norm_intensity))) %>%
  dplyr::mutate(Sample_type = factor(Sample_type, 
                                     levels = c("DSK","SK", "K")))

### Plot barplot ####
ggpubr::ggbarplot(psite_reg_count, 
                  x = "Sample_type", y = "Count", 
                  add = c("jitter", "mean_sd"),
                  label = T,
                  lab.pos = "in",
                  lab.vjust = 2.8,
                  lab.nb.digits = 0,
                  lab.size = 5,
                  fill = "Sample_type",
                  ncol = 1,
                  position = position_dodge(0.8),
                  size = 1,
                  alpha = 0.8,
                  palette = heatcolors2[[1]],
                  legend = "none",
                  x.text.angle = 45,
                  font.x = c(18, "bold"),
                  font.y = c(18, "bold"),
                  font.main = c(18, "bold"),
                  font.tickslab = c(16, "bold"),
                  title = "Class I reg. phosphosite count",
                  subtitle = "sites showing regulation in \u2265 1 construct ",
                  xlab = "Construct",
                  ylab = "Count")

ggsave(
  paste0(QC_path,
         "Barplot_validvalcount_",
         "enriched_norm_log2_intensities_L_only_regulated_woDS",
         "_120x150px_",
         projectname,
         "_",
         format(Sys.time(), '%Y%m%d_%H%M%S'),
         ".png"
  ),
  width = 150,
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
UpSetR::upset(psite_reg_upset, 
      nsets = 12,
      text.scale = 2)
dev.off()

## Create heatmap for all regulated psites ####

### Extract only intensities####

#### w/o imputation
psite_reg_ma <-psite_regulated %>%
  dplyr::select(contains("norm_intensity"))

#### with imputation
psite_reg_ma_impute <-psite_regulated %>%
  dplyr::select(contains("norm_imputed_intensity"))

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
wbook <- createWorkbook("P-Sites/Output/P-sites_ID-Count_table")
addWorksheet(wbook, "ID-counts")
writeData(wbook, sheet = 1, count_tbl)
saveWorkbook(wbook,
             paste0(DE_path, "ID-count_table_",
                    projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),
                    ".xlsx"),
             overwrite = TRUE)


###########################################################################
# Annotation-enrichment analysis###########################################
# Perform GO- and KEGG overrepresentation with the clusterProfiler 
# package
###########################################################################

## Load packages#####
library(clusterProfiler) #load clusterProfiler package
library(org.Mm.eg.db) # load annotation database for Mus musculus
library(DOSE) # load DOSE package for parsing results
library(enrichplot) # needed for GO plots

## Create output subfolder
if (dir.exists("P-Sites/Output/DE_analysis/Annotation_enrichment")) {
  
  cat("The folder already exists")
  
} else {
  
  dir.create("P-Sites/Output/DE_analysis/Annotation_enrichment",
             recursive = TRUE)
  
}

## Set relative output path####
DE_annot_path <- "P-Sites/Output/DE_analysis/Annotation_enrichment/"

## Optional - Load previous results ####
# If the upper part of the script, e.g. the DE enrichment analysis,
# was performed before we can also reload the data from the excel
# file here

# Load workbook
wbook <- openxlsx::loadWorkbook(file = choose_file())

#Re-extract data frames
psite_norm_table <- openxlsx::read.xlsx(xlsxFile = wbook, sheet = 1)
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
background_prot_ids <- psite_norm_table$Gene.names %>%
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
            showCategory = 25)
    ggsave(
      paste0(
        annot_path,
        "GO_Bubble_up_BP_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top25_",
        projectname,
        "_200x260mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 260,
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
            showCategory = 25)
    ggsave(
      paste0(
        DE_annot_path,
        "GO_Bubble_up_CC_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top25_",
        projectname,
        "_200x260mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 260,
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
            showCategory = 25)
    ggsave(
      paste0(
        DE_annot_path,
        "GO_Bubble_up_MF_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top25_",
        projectname,
        "_200x260mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 260,
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
            showCategory = 25)
    ggsave(
      paste0(
        DE_annot_path,
        "GO_Bubble_down_BP_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top25_",
        projectname,
        "_200x260mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 260,
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
            showCategory = 25)
    ggsave(
      paste0(
        DE_annot_path,
        "GO_Bubble_down_CC_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top25_",
        projectname,
        "_200x260mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 260,
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
            showCategory = 25)
    ggsave(
      paste0(
        DE_annot_path,
        "GO_Bubble_down_MF_simple_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top25_",
        projectname,
        "_200x260mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 260,
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
               paste0(DE_annot_path,"GO_overrep_", i, "_"
                      ,"pVal<0.05_qVal<0.2_", 
                      projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),
                      ".xlsx"),
               overwrite = TRUE)
  }


## KEGG-terms####
# First we generate the background lists (parameter: "universe") of all protein ids 
# followed by the individual set of differential enriched proteins for every 
# comparison (parameter: "gene")

### Over-representation analysis ####
#### Generate background list####
background_uniprot_ids <- psite_norm_table$Leading.proteins %>%
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
             showCategory = 25)
    ggsave(
      paste0(
        DE_annot_path,
        "KEGG_Bubble_up_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top25_",
        projectname,
        "_200x260mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 260,
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
            showCategory = 25)
    ggsave(
      paste0(
        DE_annot_path,
        "KEGG_Bubble_down_",
        gsub(" ", "", i),
        "_",
        "pVal<0.05_qVal<0.2_Top25_",
        projectname,
        "_200x260mm_",
        format(Sys.time(), '%Y%m%d_%H%M%S'),
        ".png"
      ),
      height = 260,
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
               paste0(DE_annot_path,"KEGG_overrep_", i, "_"
                      ,"pVal<0.05_qVal<0.2_Top30_", 
                      projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),
                      ".xlsx"),
               overwrite = TRUE)
}

################################################################################
# Venn of regulated psites in DSK-SK + DSK-K vs Dimethyl data ###########
################################################################################
# We want to do two comparisons:
# 1. DSK_SK up + DSK_K up (=microtubule-localized) vs. dimethyl up (day7/0 + 14/0) 
# 2. DSK_SK down + DSK_K down (=cytosol/nucleus localized) vs. dimethyl down (day7/0 + 14/0) 

## Extract siteIDs for each regulation ####
### Create vector of DSK vs SK upregulated sites ####
DSK_SK_up <- candidates_DSK_SK %>%
  dplyr::filter (logFC > 0.58) %>%
  unite(col = Site, 
        c(Gene.names, Amino.acid), 
        remove = T) %>% unite(col = Site,
                              c(Site,Position),
                              sep = "",
                              remove = T) %>%
  dplyr::select(c("Site")) %>% 
  separate_wider_delim(cols = everything(), 
                       delim = ";",
                       names = c("SiteID"),
                       too_many = "drop") %>% 
  filter(SiteID != "")

### Create vector of DSK vs K upregulated sites ####
DSK_K_up <- candidates_DSK_K %>%
  dplyr::filter (logFC > 0.58) %>%
  unite(col = Site, 
        c(Gene.names, Amino.acid), 
        remove = T) %>% unite(col = Site,
                              c(Site,Position),
                              sep = "",
                              remove = T) %>%
  dplyr::select(c("Site")) %>% 
  separate_wider_delim(cols = everything(), 
                       delim = ";",
                       names = c("SiteID"),
                       too_many = "drop") %>% 
  filter(SiteID != "")

### Merge DSK_SK and DSK_K up ####
microtubule_up <- rbind(DSK_SK_up, DSK_K_up) %>%
  distinct() %>%
  separate_wider_delim(cols = everything(), 
                       delim = ";",
                       names = c("Site"),
                       too_many = "drop") %>% 
  filter(Site != "")


### Create vector of DSK vs SK downregulated sites ####
DSK_SK_down <- candidates_DSK_SK %>%
  dplyr::filter (logFC < -0.58) %>%
  unite(col = Site, 
        c(Gene.names, Amino.acid), 
        remove = T) %>% unite(col = Site,
                              c(Site,Position),
                              sep = "",
                              remove = T) %>%
  dplyr::select(c("Site")) %>% 
  separate_wider_delim(cols = everything(), 
                       delim = ";",
                       names = c("SiteID"),
                       too_many = "drop") %>% 
  filter(SiteID != "")

### Create vector of DSK vs K downregulated sites ####
DSK_K_down <- candidates_DSK_K %>%
  dplyr::filter (logFC < -0.58) %>%
  unite(col = Site, 
        c(Gene.names, Amino.acid), 
        remove = T) %>% unite(col = Site,
                              c(Site,Position),
                              sep = "",
                              remove = T) %>%
  dplyr::select(c("Site")) %>% 
  separate_wider_delim(cols = everything(), 
                       delim = ";",
                       names = c("SiteID"),
                       too_many = "drop") %>% 
  filter(SiteID != "")

### Merge DSK_SK and DSK_K down ####
cytonuc_up <- rbind(DSK_SK_down, DSK_K_down) %>%
  distinct() %>%
  separate_wider_delim(cols = everything(), 
                       delim = ";",
                       names = c("Site"),
                       too_many = "drop") %>% 
  filter(Site != "")


### dimethyl up = psites upregulated in day7/0 and/or day14/0 ####
DML_7_0_up <- DML_7_0_all %>% 
  dplyr::filter (`2Xup` == "+") %>%
  dplyr::select(c("Gene.names", "Amino.acid", "Position")) %>%
  unite(col = Site, 
        c(Gene.names, Amino.acid), 
        remove = T) %>% unite(col = Site,
                              c(Site,Position),
                              sep = "",
                              remove = T) 

DML_14_0_up <- DML_14_0_all %>% 
  dplyr::filter (`2Xup` == "+") %>%
  dplyr::select(c("Gene.names", "Amino.acid", "Position")) %>%
  unite(col = Site, 
        c(Gene.names, Amino.acid), 
        remove = T) %>% unite(col = Site,
                              c(Site,Position),
                              sep = "",
                              remove = T) 

DML_up <- rbind(DML_7_0_up,DML_14_0_up) %>%
  distinct() %>%
  separate_wider_delim(cols = everything(), 
                       delim = ";",
                       names = c("Site"),
                       too_many = "drop") %>% 
  filter(Site != "")

### dimethyl down = psites downregulated in day7/0 and/or day14/0 ####
DML_7_0_down <- DML_7_0_all %>% 
  dplyr::filter (`2Xdown` == "+") %>%
  dplyr::select(c("Gene.names", "Amino.acid", "Position")) %>%
  unite(col = Site, 
        c(Gene.names, Amino.acid), 
        remove = T) %>% unite(col = Site,
                              c(Site,Position),
                              sep = "",
                              remove = T) 

DML_14_0_down <- DML_14_0_all %>% 
  dplyr::filter (`2Xdown` == "+") %>%
  dplyr::select(c("Gene.names", "Amino.acid", "Position")) %>%
  unite(col = Site, 
        c(Gene.names, Amino.acid), 
        remove = T) %>% unite(col = Site,
                              c(Site,Position),
                              sep = "",
                              remove = T) 

DML_down <- rbind(DML_7_0_down,DML_14_0_down) %>%
  distinct() %>%
  separate_wider_delim(cols = everything(), 
                       delim = ";",
                       names = c("Site"),
                       too_many = "drop") %>% 
  filter(Site != "")

## Create lists for each Venn to be ploted ####
venn_a <- list(microtubule_up = microtubule_up$Site,
               DML_up = DML_up$Site)

venn_b <- list(cytonuc_up = cytonuc_up$Site,
               DML_down = DML_down$Site)

## Calculate Venn overlaps & write to file
venn_a_counts <- get.venn.partitions(venn_a)
write.xlsx(x = venn_a_counts, 
           file = paste0("P-Sites/Output/Venn_diagrams/Venn-table_microtubules-up_vs_DML-up_",
                         projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),".xlsx"))

venn_b_counts <- get.venn.partitions(venn_b)
write.xlsx(x = venn_b_counts, 
           file = paste0("P-Sites/Output/Venn_diagrams/Venn-table_cytosol-nucleus-up_vs_DML-down_",
                         projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),".xlsx"))

## Plot Venn diagrams ####
venn.diagram(venn_a,
             filename = paste0("P-Sites/Output/Venn_diagrams/Venn_microtubules-up_vs_DML-up_",
                               projectname, "_",format(Sys.time(), 
                                                       '%Y%m%d_%H%M%S'),".png"),
             category.names = c(microtubule_up = "Microtubules up",
                                DML_up = "Dimethyl up"),
             print.mode = c("raw", "percent"),
             sigdigs = 1,
             fill = c("#EFC000FF", "#5CB85CFF"),
             fontfamily = "sans",
             cat.fontfamily = "sans",
             cat.fontface = "bold",
             cat.col = c("#EFC000FF", "#5CB85CFF"),
             cat.cex = 1.4,
             imagetype = "png",
             cat.just = list(a = c(-0.5, -4),
                             b = c(1.7, -4)),
             margin = 0.1,
             euler = F,
             scale = F)

venn.diagram(venn_b,
             filename = paste0("P-Sites/Output/Venn_diagrams/Venn_cytosol-nucleus-down_vs_DML-down__",
                               projectname, "_",format(Sys.time(), 
                                                       '%Y%m%d_%H%M%S'),".png"),
             category.names = c(cytonuc_up = "Cytosol/nucleus up",
                                DML_down = "Dimethyl down"),
             print.mode = c("raw", "percent"),
             sigdigs = 1,
             fill = c("#D43F3AFF", "#9632B8FF"),
             fontfamily = "sans",
             cat.fontfamily = "sans",
             cat.fontface = "bold",
             cat.col = c("#D43F3AFF", "#9632B8FF"),
             cat.cex = 1.4,
             imagetype = "png",
             cat.just = list(a = c(-0.5, -4),
                             b = c(1.7, -4)),
             margin = 0.1,
             euler = F,
             scale = F)
################################################################################
# Write session info to file ###################################################
################################################################################

#Save session_info to file for traceabillity & transparency reasons####
devtools::session_info(to_file = paste0("P-Sites/Session_info_",
                                        projectname, "_", 
                                        format(Sys.time(),'%Y%m%d_%H%M%S'),
                                        ".txt"))