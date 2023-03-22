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
    "ggplot2",
    "dplyr",
    "stringr",
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
    "RColorBrewer"
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

BiocManager::install(c('limma', "clusterProfiler",
                       "MsCoreUtils", "pcaMethods", "org.Mm.eg.db"))

library("MsCoreUtils")
library(clusterProfiler) #load clusterProfiler package
library(org.Mm.eg.db) # load annotation database for Mus musculus
library(DOSE) # load DOSE package for parsing results
library(enrichplot) # needed for GO plots

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

# Perform GO-overrepresentation and GSEA analyses with the clusterProfiler 
# package


## Create output subfolder
if (dir.exists("pooled_annotation_enrichment")) {
  
  cat("The folder already exists")
  
} else {
  
  dir.create("pooled_annotation_enrichment")
  
}

## Set relative output path####
annot_path <- "pooled_annotation_enrichment/"

## Set project name ####
# This name is added to all result tables and plots. 
#Please change for every project/dataset to be analyzed.
projectname <- "03_BioID-DCLK_proteome"

##Load previous results ####
# If the  DE enrichment analysis,
# was performed before, we can also reload the data from the excel
# file here

# Load workbook
wbook <- openxlsx::loadWorkbook(file = choose_file())

#Re-extract data frames
prot_table <- openxlsx::read.xlsx(xlsxFile = wbook, sheet = 1)
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
background_prot_ids <- prot_table$Gene.names %>%
  str_split(";", simplify = TRUE) %>%
  as.data.frame()  %>% 
  dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
background_prot_ids <- background_prot_ids$V1
background_prot_ids <- na.omit(background_prot_ids)

#### Generate DE-protein lists####

# DSK vs SK
DSK_SK_up <- filter(candidates_DSK_SK,
                    logFC > 0.58 & adj.P.Val < 0.05 ) %>%
  dplyr::select(c("Gene.names"))
DSK_SK_up <- str_split(DSK_SK_up$Gene.names, ";", simplify = TRUE) %>%
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
DSK_SK_up <- DSK_SK_up$V1
DSK_SK_up <- na.omit(DSK_SK_up)
  
DSK_SK_down <- filter(candidates_DSK_SK,
                      logFC < -0.58 & adj.P.Val < 0.05 ) %>%
  dplyr::select(c("Gene.names"))
DSK_SK_down <- str_split(DSK_SK_down$Gene.names, ";", simplify = TRUE) %>%
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
DSK_SK_down <- DSK_SK_down$V1
DSK_SK_down <- na.omit(DSK_SK_down)

#DSK vs K
DSK_K_up <- filter(candidates_DSK_K,
                    logFC > 0.58 & adj.P.Val < 0.05 ) %>%
  dplyr::select(c("Gene.names"))
DSK_K_up <- str_split(DSK_K_up$Gene.names, ";", simplify = TRUE) %>%
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
DSK_K_up <- DSK_K_up$V1
DSK_K_up <- na.omit(DSK_K_up)

DSK_K_down <- filter(candidates_DSK_K,
                      logFC < -0.58 & adj.P.Val < 0.05 ) %>%
  dplyr::select(c("Gene.names"))
DSK_K_down <- str_split(DSK_K_down$Gene.names, ";", simplify = TRUE) %>%
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
DSK_K_down <- DSK_K_down$V1
DSK_K_down <- na.omit(DSK_K_down)

#### Combine DSK_SK & DSK_K ####
DSK_SK_K_up <- unique(c(DSK_SK_up,DSK_K_up))
DSK_SK_K_down <- unique(c(DSK_SK_down,DSK_K_down))


#### Perform over-representation analyses####
  ego_up_BP <- enrichGO(
    gene          = DSK_SK_K_up,
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
    gene          = DSK_SK_K_up,
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
    gene          = DSK_SK_K_up,
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
    gene          = DSK_SK_K_down,
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
    gene          = DSK_SK_K_down,
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
    gene          = DSK_SK_K_down,
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
  
  if(nrow(base::as.data.frame(ego_up_BP)) > 0)
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
        "DSK_SK_K",
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
        "DSK_SK_K",
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
        "DSK_SK_K",
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
        "DSK_SK_K",
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
        "DSK_SK_K",
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
        "DSK_SK_K",
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
               paste0(annot_path,"GO_overrep_", "DSK_SK_K", "_"
                      ,"pVal<0.05_qVal<0.2_Top30_", 
                      projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),
                      ".xlsx"),
               overwrite = TRUE)


## KEGG-terms####
# First we generate the background lists (parameter: "universe") of all protein ids 
# followed by the individual set of differential enriched proteins for every 
# comparison (parameter: "gene")

### Over-representation analysis ####
#### Generate background list####
background_uniprot_ids <- prot_table$Majority.protein.IDs %>%
  str_split(";", simplify = TRUE) %>%
  as.data.frame()  %>% 
  dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
background_uniprot_ids <- background_uniprot_ids$V1
background_uniprot_ids <- na.omit(background_uniprot_ids)

#### Generate DE-protein lists####
# DSK vs SK
DSK_SK_up_uniprot_ids <- filter(candidates_DSK_SK,
                    logFC > 0.58 & adj.P.Val < 0.05 ) %>%
  dplyr::select(c("Majority.protein.IDs"))
DSK_SK_up_uniprot_ids <- str_split(DSK_SK_up_uniprot_ids$Majority.protein.IDs, ";", simplify = TRUE) %>%
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
DSK_SK_up_uniprot_ids <- DSK_SK_up_uniprot_ids$V1
DSK_SK_up_uniprot_ids <- na.omit(DSK_SK_up_uniprot_ids)

DSK_SK_down_uniprot_ids <- filter(candidates_DSK_SK,
                      logFC < -0.58 & adj.P.Val < 0.05 ) %>%
  dplyr::select(c("Majority.protein.IDs"))
DSK_SK_down_uniprot_ids <- str_split(DSK_SK_down_uniprot_ids$Majority.protein.IDs, ";", simplify = TRUE) %>%
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
DSK_SK_down_uniprot_ids <- DSK_SK_down_uniprot_ids$V1
DSK_SK_down_uniprot_ids <- na.omit(DSK_SK_down_uniprot_ids)

#DSK vs K
DSK_K_up_uniprot_ids <- filter(candidates_DSK_K,
                   logFC > 0.58 & adj.P.Val < 0.05 ) %>%
  dplyr::select(c("Majority.protein.IDs"))
DSK_K_up_uniprot_ids <- str_split(DSK_K_up_uniprot_ids$Majority.protein.IDs, ";", simplify = TRUE) %>%
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
DSK_K_up_uniprot_ids <- DSK_K_up_uniprot_ids$V1
DSK_K_up_uniprot_ids <- na.omit(DSK_K_up_uniprot_ids)

DSK_K_down_uniprot_ids <- filter(candidates_DSK_K,
                     logFC < -0.58 & adj.P.Val < 0.05 ) %>%
  dplyr::select(c("Majority.protein.IDs"))
DSK_K_down_uniprot_ids <- str_split(DSK_K_down_uniprot_ids$Majority.protein.IDs, ";", simplify = TRUE) %>%
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
DSK_K_down_uniprot_ids <- DSK_K_down_uniprot_ids$V1
DSK_K_down_uniprot_ids <- na.omit(DSK_K_down_uniprot_ids)

#### Combine DSK_SK & DSK_K ####
DSK_SK_K_up_uniprot_ids <- unique(c(DSK_SK_up_uniprot_ids,DSK_K_up_uniprot_ids))
DSK_SK_K_down_uniprot_ids <- unique(c(DSK_SK_down_uniprot_ids,DSK_K_down_uniprot_ids))


  #### Perform over-representation analyses####
  
  ego_down_KEGG <- enrichKEGG(
    gene          = DSK_SK_K_down_uniprot_ids,
    universe      = background_uniprot_ids,
    keyType = "uniprot",
    organism       = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff  = 0.2
  )
  
  ego_up_KEGG <- enrichKEGG(
    gene          = DSK_SK_K_up_uniprot_ids,
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
        "DSK_SK_K",
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
        "DSK_SK_K",
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
               paste0(annot_path,"KEGG_overrep_", "DSK_SK_K", "_"
                      ,"pVal<0.05_qVal<0.2_Top30_", 
                      projectname,"_",format(Sys.time(), '%Y%m%d_%H%M%S'),
                      ".xlsx"),
               overwrite = TRUE)
