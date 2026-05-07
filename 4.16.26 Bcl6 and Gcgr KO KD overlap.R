# Initialize --------------------------------------------------------------

## Install + load required packages

# install.packages("BiocManager")
# install.packages("tidyverse")
# install.packages("ggVennDiagram")

library("pacman")

pacman::p_load(
  BiocManager,
  GO.db, 
  ChIPseeker, 
  TxDb.Mmusculus.UCSC.mm10.knownGene, 
  org.Mm.eg.db,
  clusterProfiler,
  DESeq2,
  biomaRt,
  enrichplot,
  tidyverse,
  ggVennDiagram,
  dplyr,
  readxl,
  R.utils,
  EnhancedVolcano,
  try.bioconductor = TRUE
  )

# ChIP Section ------------------------------------------------------------

# Load bcl6 peak file - using fed
bcl6_fed <- readPeakFile("GSE118788_BCL6-C57-Fed-peaks.bed.txt")
bcl6_fast <- readPeakFile("GSE118788_BCL6-C57-Fast-peaks.bed.txt")

# Processing function
process_bcl6 <- function(peak_file) {
  # Annotate peaks
  anno <- annotatePeak(
    peak_file,
    TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
    tssRegion = c(-3000, 3000),
    annoDb = "org.Mm.eg.db"
  )
  # Extract gene symbols
  ids <- as.data.frame(anno)$geneId

  symbols <- bitr(
    ids,
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = org.Mm.eg.db
  )$SYMBOL
  
  targets <- as.data.frame(anno) %>%
    pull(ENSEMBL) %>%
    unique()

  return(list(anno = anno, genes = ids, symbols = symbols, targets = targets))
}

# Run for FAST
fast_results <- process_bcl6(bcl6_fast)
peak_anno_fast <- fast_results$anno
bcl6_genes_fast <- fast_results$genes
bcl6_symbols_fast <- fast_results$symbols
bcl6_targets_fast <- fast_results$targets

# Run for FED
fed_results <- process_bcl6(bcl6_fed)
peak_anno_fed <- fed_results$anno
bcl6_genes_fed <- fed_results$genes
bcl6_symbols_fed <- fed_results$symbols
bcl6_targets_fed <- fed_results$targets


# RNASeq Section - old ----------------------------------------------------------

# # 1: Liver-specific gcgr KO vs WT (counts), M+F, Fast unknown
# counts <- read.delim(
#   gzfile("GSE302496_counts.txt.gz"),
#   comment.char = "#",
#   header = TRUE,
#   sep = "\t",
#   check.names = FALSE
# )
# 
# rownames(counts) <- counts$Geneid
# counts <- counts[, !(colnames(counts) %in% c("Geneid", "Chr", "Start", "End", "Strand", "Length"))]
# 
# #2: Al Powers gcgr blocking mAb 10 days (relative expression), M+F, Fasted
# gunzip("GSE89035_normalized_expression_values_all_2552ACP.xls.gz", overwrite = FALSE)
# expr <- read_excel("GSE89035_normalized_expression_values_all_2552ACP.xls")
# 
# str(expr)
# head(expr)
# rownames(expr) <- expr[[1]]
# expr <- expr[,-1]
# 
# 
# mart <- useEnsembl(
#   biomart = "genes",
#   dataset = "mmusculus_gene_ensembl"
# )
# 
# map <- getBM(
#   attributes = c("ensembl_gene_id", "mgi_symbol"),
#   mart = mart
# )
# 
# de_genes <- expr$symbol[expr$padj < 0.05 & abs(expr$log2FoldChange) > 1]


# RNASeq Section - new ----------------------------------------------------

# deg_data <- read_excel("NIHMS880297-supplement-3.xlsx", sheet = "Significantly altered in both", skip = 1)
# my_degs <- deg_data %>%
#   pull(`Ensembl ID`) %>%
#   unique()

gcgr_KO_data <- read_excel("NIHMS880297-supplement-3.xlsx", sheet = "FC-Gcgr- vs Gcgr+", skip = 1)
gcgr_mAbKD_data <- read_excel("NIHMS880297-supplement-3.xlsx", sheet = "FC-Gcgr mAb vs Gcgr+", skip = 1)

gcgr_KO_degs <- gcgr_KO_data %>%
  pull(`Ensembl ID`) %>%
  unique()

gcgr_mAbKD_degs <- gcgr_mAbKD_data %>%
  pull(`Ensembl ID`) %>%
  unique()

# Identify Overlap -----------------------------------------------------------------

# overlap_genes <- intersect(my_degs, bcl6_targets_fast)

gcgr_KO_overlap_genes <- intersect(gcgr_KO_degs, bcl6_targets_fed)
gcgr_mAbKD_overlap_genes <- intersect(gcgr_mAbKD_degs, bcl6_targets_fed)

#convert all to symbols/entrez function

get_symbols_entrez <- function(genelist) {
  
  symbols <- bitr(
    genelist,
    fromType = "ENSEMBL",
    toType = "SYMBOL",
    OrgDb = org.Mm.eg.db
  )$SYMBOL
  
  entrez <- bitr(
    genelist,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )$ENTREZID
  
  return(list(symbols = symbols, entrez = entrez))
}

gcgr_KO_results <- get_symbols_entrez(gcgr_KO_degs)
gcgr_KO_symbols <- gcgr_KO_results$symbols
gcgr_KO_entrez <- gcgr_KO_results$entrez

gcgr_mAbKD_results <- get_symbols_entrez(gcgr_mAbKD_degs)
gcgr_mAbKD_symbols <- gcgr_mAbKD_results$symbols
gcgr_mAbKD_entrez <- gcgr_mAbKD_results$entrez




# Venn Diagram ---------------------------------------------------------------

gene_lists <- list(GCGR_KO = gcgr_KO_symbols,
                   GCGR_mAbKD = gcgr_mAbKD_symbols,
                   BCL6_Modulated = bcl6_symbols_fed)

ggVennDiagram(gene_lists) + 
  scale_fill_gradient(low="blue", high = "red") +
  labs(title = "Intersection of BCL6 Binding and DEGs in GCGR KO + KD")

# Volcano Plot ------------------------------------------------------------

# add column to DEGs correlating to BCL6 peaks
gcgr_KO_data <- gcgr_KO_data %>%
  mutate(is_bcl6_target = if_else(`Ensembl ID` %in% bcl6_targets_fed, "Y", "N"))
gcgr_mAbKD_data <- gcgr_mAbKD_data %>%
  mutate(is_bcl6_target = if_else(`Ensembl ID` %in% bcl6_targets_fed, "Y", "N"))


# rename columns in simpler/easier to use ways

gcgr_KO_data <- gcgr_KO_data %>% rename(
                                        "symbol" = `Gene Symbol`,
                                        "ensembl_id" = `Gene ID`,
                                        "FC_KO" = `FC-Gcgr-/- vs Gcgr+/+`,
                                        "Log_FC_KO" = `Log FC-Gcgr-/- vs Gcgr+/+`,
                                        "ur_dr_KO" = `Regulation FC-Gcgr-/- vs Gcgr+/+`,
                                        "pval" = `p-value ( Z Test )`,
                                        "pval_adj" = `Corrected p-value ( Z Test )`,
                                        "entrez_id" = `Entrez ID`
                                        )


gcgr_mAbKD_data <- gcgr_mAbKD_data %>% rename(
                                              "symbol" = `Gene Symbol`,
                                              "ensembl_id" = `Gene ID`,
                                              "FC_KD" = `FC-FC-Gcgr mAb vs Gcgr+/+`,
                                              "Log_FC_KD" = `Log FC-Gcgr mAb vs Gcgr+/+`,
                                              "ur_dr_KD" = `Regulation FC-Gcgr mAb vs Gcgr+/+`,
                                              "pval" = `p-value ( Z Test )`,
                                              "pval_adj" = `Corrected p-value ( Z Test )`,
                                              "entrez_id" = `Entrez ID`
                                              )


# remove anything that isn't modulated by bcl6, and create proper rownames for volcano plot
KO_data_clean <- gcgr_KO_data %>%
  filter(is_bcl6_target == "Y") %>%
  column_to_rownames(var = "symbol")
KD_data_clean <- gcgr_mAbKD_data %>%
  filter(is_bcl6_target == "Y") %>%
  column_to_rownames(var = "symbol")


#from scratch plot function
volcano_plotter <- function(df, pval_limit, logFC_limit, title) {
  
  pval_cutoff <- -log10(pval_limit)
  logFC_cutoff <- logFC_limit
  
  df <- df %>%
    mutate(
      plot_pval_adj = ifelse(pval_adj == 0, 1e-17, pval_adj),
      plot_pval = ifelse(pval == 0, 1e-17, pval),
      
      label = ifelse(
        rank(Log_FC, ties.method = "first") <= 15 |
          rank(-Log_FC, ties.method = "first") <= 15,
        rownames(df),
        NA
      ),
      
      plot_group = ifelse(
        abs(Log_FC) < logFC_cutoff |
          -log10(plot_pval_adj) < pval_cutoff,
        "ns",
        ur_dr
      )
    )
  
  ggplot(df, aes(x = Log_FC, y = -log10(plot_pval_adj), col = plot_group)) +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), col = "gray", linetype = "dashed") +
    geom_hline(yintercept = pval_cutoff, col = "gray", linetype = "dashed") +
    geom_point(size = 2) +
    scale_color_manual(values = c(
      "ns" = "gray",
      down = "#bb0c00",
      up = "#00AFBB"
    )) +
    labs(
      x = expression("log"[2]*"FC"),
      y = expression("-log"[10]*"p-value")
    ) +
    coord_cartesian(ylim = c(0, 20), xlim = c(-7, 7)) +
    scale_x_continuous(breaks = seq(-10, 10, 2)) +
    scale_y_continuous(breaks = seq(0, 20, 5)) +
    ggtitle(title) +
    geom_text_repel(
      aes(label = label),
      max.overlaps = Inf,
      point.padding = 0.1,
      box.padding = 0.5,
      force = 2,
      segment.color = "grey30"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
}
volcano_plotter(KO_data_clean, 0.05, 0.6, "DEG in GCGR-/- KO + BCL6 peak")
volcano_plotter(KD_data_clean, 0.05, 0.6, "DEG in GCGR mAB KD + BCL6 peak")




#EnrichGO
ego_bp <- enrichGO(
  gene          = overlap_entrez,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

head(ego_bp)

dotplot(ego_bp, showCategory = 15)
cnetplot(ego_bp, showCategory = 5)


######
#Next steps:
#Separate DEGs by KO vs KD
#Repeat analysis
#Make a volcano plot of DEGs with color if modified by BCL6
#Figure out upregulation vs downregulation and log fold effect sizes
