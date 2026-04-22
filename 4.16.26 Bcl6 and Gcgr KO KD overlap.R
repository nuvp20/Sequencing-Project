# Initialize --------------------------------------------------------------

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
  usethis,
  gitcreds,
  data.table,
  rtracklayer,
  try.bioconductor = TRUE
  )


# Bcl6 ChIP Section ------------------------------------------------------------

# Load bcl6 peak file - using fed
bcl6_fed <- readPeakFile("GSE118788_BCL6-C57-Fed-peaks.bed.txt")
bcl6_fast <- readPeakFile("GSE118788_BCL6-C57-Fast-peaks.bed.txt")

# Processing function
process_chip <- function(peak_file) {
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
fast_results <- process_chip(bcl6_fast)
peak_anno_fast <- fast_results$anno
bcl6_genes_fast <- fast_results$genes
bcl6_symbols_fast <- fast_results$symbols
bcl6_targets_fast <- fast_results$targets

# Run for FED
fed_results <- process_chip(bcl6_fed)
peak_anno_fed <- fed_results$anno
bcl6_genes_fed <- fed_results$genes
bcl6_symbols_fed <- fed_results$symbols
bcl6_targets_fed <- fed_results$targets


# RNASeq Section - new ----------------------------------------------------

# Read files with gcgr info
gcgr_KO_data <- read_excel("NIHMS880297-supplement-3.xlsx", sheet = "FC-Gcgr- vs Gcgr+", skip = 1)
gcgr_mAbKD_data <- read_excel("NIHMS880297-supplement-3.xlsx", sheet = "FC-Gcgr mAb vs Gcgr+", skip = 1)

# Pull each unique DEG
gcgr_KO_degs <- gcgr_KO_data %>%
  pull(`Ensembl ID`) %>%
  unique()

gcgr_mAbKD_degs <- gcgr_mAbKD_data %>%
  pull(`Ensembl ID`) %>%
  unique()

# Identify Overlap -----------------------------------------------------------------

# Find the overlap between Bcl6 and each Gcgr condition
gcgr_KO_overlap_genes <- intersect(gcgr_KO_degs, bcl6_targets_fed)
gcgr_mAbKD_overlap_genes <- intersect(gcgr_mAbKD_degs, bcl6_targets_fed)

# Function - given ensembl, pull symbol and entrez ID
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

#Pull alternate IDs for all DEGs
gcgr_KO_results <- get_symbols_entrez(gcgr_KO_degs)
gcgr_KO_symbols <- gcgr_KO_results$symbols
gcgr_KO_entrez <- gcgr_KO_results$entrez

gcgr_mAbKD_results <- get_symbols_entrez(gcgr_mAbKD_degs)
gcgr_mAbKD_symbols <- gcgr_mAbKD_results$symbols
gcgr_mAbKD_entrez <- gcgr_mAbKD_results$entrez


# Venn Diagram ---------------------------------------------------------------

# All DE/Modulated genes, visualize magnitude of overlap
gene_lists <- list(GCGR_KO = gcgr_KO_symbols,
                   GCGR_mAbKD = gcgr_mAbKD_symbols,
                   BCL6_Modulated = bcl6_symbols_fed)

ggVennDiagram(gene_lists) + 
  scale_fill_gradient(low="blue", high = "red") +
  labs(title = "Intersection of Bcl6 Binding and DEGs in Gcgr KO + KD") +
  scale_x_continuous(expand = expansion(mult = 0.2))

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
                                        "FC" = `FC-Gcgr-/- vs Gcgr+/+`,
                                        "Log_FC" = `Log FC-Gcgr-/- vs Gcgr+/+`,
                                        "ur_dr" = `Regulation FC-Gcgr-/- vs Gcgr+/+`,
                                        "pval" = `p-value ( Z Test )`,
                                        "pval_adj" = `Corrected p-value ( Z Test )`,
                                        "entrez_id" = `Entrez ID`
                                        )

gcgr_mAbKD_data <- gcgr_mAbKD_data %>% rename(
                                              "symbol" = `Gene Symbol`,
                                              "ensembl_id" = `Gene ID`,
                                              "FC" = `FC-FC-Gcgr mAb vs Gcgr+/+`,
                                              "Log_FC" = `Log FC-Gcgr mAb vs Gcgr+/+`,
                                              "ur_dr" = `Regulation FC-Gcgr mAb vs Gcgr+/+`,
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


# Over-representation analysis --------------------------------------------

ego_KO <- enrichGO(
  gene = gcgr_KO_overlap_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "ALL", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

dotplot(ego_KO,
        title = "Pathways enriched in overlap of Bcl6 + Gcgr -/- KO"
        )

ego_KD <- enrichGO(
  gene = gcgr_mAbKD_overlap_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "ALL", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

dotplot(ego_KD, 
        title = "Pathways enriched in overlap of Bcl6 + Gcgr mAb KD"
        )




# Creb ChIP Section -------------------------------------------------------

creb_peaks_fast <- import("GSE45674_CREB_peaks_fasted_final.bed.gz")
creb_fast_results <- process_chip(creb_peaks_fast)
creb_peak_anno_fast <- creb_fast_results$anno
creb_genes_fast <- creb_fast_results$genes
creb_symbols_fast <- creb_fast_results$symbols
creb_targets_fast <- creb_fast_results$targets

creb_peaks_fed <- import("GSE45674_CREB_peaks_refed_final.bed.gz")
creb_fed_results <- process_chip(creb_peaks_fed)
creb_peak_anno_fed <- creb_fed_results$anno
creb_genes_fed <- creb_fed_results$genes
creb_symbols_fed <- creb_fed_results$symbols
creb_targets_fed <- creb_fed_results$targets

# Bcl6-Creb overlap + Venn diagram -------------------------------------------------------

bcl6_creb_overlap_fast <- intersect(bcl6_targets_fast, creb_targets_fast)
bcl6_creb_overlap_fed <- intersect(bcl6_targets_fed, creb_targets_fed)

gene_lists <- list(Bcl6_fed <- bcl6_targets_fed,
                   Bcl6_fasted <- bcl6_targets_fast,
                   Creb_fed <- creb_targets_fed,
                   Creb_fasted <- creb_targets_fast)

ggVennDiagram(gene_lists, category.names = c("Bcl6_fed", "Bcl6_fasted", "Creb_fed", "Creb_fasted")) + 
  scale_fill_gradient(low="blue", high = "red") +
  labs(title = "Intersection of Bcl6 and Creb") +
  scale_x_continuous(expand = expansion(mult = 0.2))

