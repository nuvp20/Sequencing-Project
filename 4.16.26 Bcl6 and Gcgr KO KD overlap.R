# Initialize --------------------------------------------------------------
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
  patchwork,
  gridGraphics,
  gridExtra,
  ggpubr,
  eulerr,
  janitor,
  openxlsx,
  ReactomePA,
  KEGGREST,
  FELLA,
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

kdkofit <- euler(gene_lists)

kdkovenn <- plot(
                  kdkofit,
                  quantities = TRUE,
                  fills = list(alpha = 0.6),
                  edges = TRUE,
                  labels = list(font = 2),
                  main = "Intersection of Bcl6 Binding and DEGs in Gcgr KO + KD"
                )
kdkovenn

# Not using ggVennDiagram
# ggVennDiagram(gene_lists) + 
#   scale_fill_gradient(low="blue", high = "red") +
#   labs(title = "Intersection of Bcl6 Binding and DEGs in Gcgr KO + KD") +
#   scale_x_continuous(expand = expansion(mult = 0.2))

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
volcano_plotter <- function(df, pval_limit, logFC_limit, ymax = 20, xlims = 7, title, rownames = TRUE) {
  
  if(rownames) {
    df <- df %>% rownames_to_column(var = "symbol")
  }
  
  pval_cutoff <- -log10(pval_limit)
  logFC_cutoff <- logFC_limit
  
  df <- df %>%
    mutate(
      plot_pval_adj = ifelse(pval_adj == 0, 1e-17, pval_adj),
      plot_pval = ifelse(pval == 0, 1e-17, pval),
      
      label = ifelse(
        rank(Log_FC, ties.method = "first") <= 15 |
          rank(-Log_FC, ties.method = "first") <= 15,
        symbol,
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
    coord_cartesian(ylim = c(0, ymax), xlim = c(-xlims, xlims)) +
    scale_x_continuous(breaks = seq(-(xlims + 2), (xlims + 2), 2)) +
    scale_y_continuous(breaks = seq(0, ymax, ymax/20)) +
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

# some values are missing????
vp1 <- volcano_plotter(KO_data_clean, 0.05, 0.6, title = "DEG in GCGR-/- KO + BCL6 peak")
vp2 <- volcano_plotter(KD_data_clean, 0.05, 0.6, title = "DEG in GCGR mAB KD + BCL6 peak")
vp1
vp2
volcanoes <- vp1 + vp2
volcanoes

# Over-representation analysis --------------------------------------------

ego_KO1 <- enrichGO(
  gene = gcgr_KO_overlap_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
ego_KO2 <- enrichGO(
  gene = gcgr_KO_overlap_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "MF", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
ego_KO3 <- enrichGO(
  gene = gcgr_KO_overlap_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "CC", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

kogo1 <- dotplot(ego_KO1, title = "Bcl6 + Gcgr -/- KO (Biological Process)") + scale_y_discrete(labels = function(x) str_wrap(x, width = 60))
kogo2 <- dotplot(ego_KO2, title = "Bcl6 + Gcgr -/- KO (Molecular Function)") + scale_y_discrete(labels = function(x) str_wrap(x, width = 60))
kogo3 <- dotplot(ego_KO3, title = "Bcl6 + Gcgr -/- KO (Cellular Component)") + scale_y_discrete(labels = function(x) str_wrap(x, width = 60))
kogoplot <- kogo1 / kogo2 / kogo3 + plot_layout(guides = "collect") & theme(legend.position = "right", axis.text.y = element_text(size = 8))
kogoplot

ego_KD1 <- enrichGO(
  gene = gcgr_mAbKD_overlap_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
ego_KD2 <- enrichGO(
  gene = gcgr_mAbKD_overlap_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "MF", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
ego_KD3 <- enrichGO(
  gene = gcgr_mAbKD_overlap_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "CC", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

kdgo1 <- dotplot(ego_KD1, title = "Bcl6 + Gcgr mAb KD (Biological Process)") + scale_y_discrete(labels = function(x) str_wrap(x, width = 60))
kdgo2 <- dotplot(ego_KD2, title = "Bcl6 + Gcgr mAb KD (Molecular Function)") + scale_y_discrete(labels = function(x) str_wrap(x, width = 60))

# No CC enrichment, so kdgo3 gets special handling
if (nrow(as.data.frame(ego_KD3)) == 0) {
  kdgo3 <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No significant CC enrichment",
             size = 5) +
    theme_void() +
    ggtitle("Bcl6 + Gcgr mAb KD (Cellular Component)")
} else {
  kdgo3 <- dotplot(ego_KD3, showCategory = 10,
                   title = "Bcl6 + Gcgr mAb KD (Cellular Component)") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40))
}

kdgoplot <- kdgo1 / kdgo2 / kdgo3 + plot_layout(guides = "collect") & theme(legend.position = "right", axis.text.y = element_text(size = 8))
kdgoplot

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

# Bcl6-Creb overlap -------------------------------------------------------

bcl6_creb_overlap_fast <- intersect(bcl6_targets_fast, creb_targets_fast)
bcl6_creb_overlap_fed <- intersect(bcl6_targets_fed, creb_targets_fed)

# Creb Venn Diagram (OLD) -------------------------------------------------------

# gene_lists <- list(Bcl6_fed <- bcl6_targets_fed,
#                    Bcl6_fasted <- bcl6_targets_fast,
#                    Creb_fed <- creb_targets_fed,
#                    Creb_fasted <- creb_targets_fast)
# 
# ggVennDiagram(gene_lists, category.names = c("Bcl6_fed", "Bcl6_fasted", "Creb_fed", "Creb_fasted")) + 
#   scale_fill_gradient(low="blue", high = "red") +
#   labs(title = "Intersection of Bcl6 and Creb") +
#   scale_x_continuous(expand = expansion(mult = 0.2))
# 
# 
# fed_lists <- list (Bcl6_fed <- bcl6_targets_fed,
#                    Creb_fed <- creb_targets_fed)
# 
# fed <- ggVennDiagram(fed_lists, category.names = c("Bcl6 fed", "Creb fed")) +
#   scale_fill_gradient(low = "blue", high = "red") + 
#   labs(title = "Bcl6 + Creb overlap, fed") + 
#   scale_x_continuous(expand = expansion(mult = 0.2)) +
#   coord_flip()
# fed
# 
# fast_lists <- list (Bcl6_fast <- bcl6_targets_fast,
#                    Creb_fast <- creb_targets_fast)
# 
# fast <- ggVennDiagram(fast_lists, category.names = c("Bcl6 fasted", "Creb fasted")) +
#   scale_fill_gradient(low = "blue", high = "red") + 
#   labs(title = "Bcl6 + Creb overlap, fasted") + 
#   scale_x_continuous(expand = expansion(mult = 0.2)) +
#   coord_flip()
# fast

# Creb Venn Diagram (NEW) -------------------------------------------------

gene_lists <- list(
  Bcl6_fed = bcl6_targets_fed,
  Bcl6_fasted = bcl6_targets_fast,
  Creb_fed = creb_targets_fed,
  Creb_fasted = creb_targets_fast
)

fed_lists <- list(
  Bcl6_fed = bcl6_targets_fed,
  Creb_fed = creb_targets_fed
)

fast_lists <- list(
  Bcl6_fasted = bcl6_targets_fast,
  Creb_fasted = creb_targets_fast
)

fit_all <- euler(gene_lists)
fit_fed <- euler(fed_lists)
fit_fast <- euler(fast_lists)

creb_all_venn <- plot(
                      fit_all,
                      quantities = TRUE,
                      fills = list(alpha = 0.6),
                      edges = TRUE,
                      labels = list(font = 2),
                      main = "Intersection of Bcl6 and Creb (Fed and Fasted)"
                    )

creb_fed_venn <- plot(
                      fit_fed,
                      quantities = TRUE,
                      fills = list(alpha = 0.6, col = c("#4C78A8", "#F58518")),
                      edges = TRUE,
                      labels = list(font = 2),
                      main = "Bcl6 + Creb overlap (Fed)"
                    )

creb_fast_venn <- plot(
                        fit_fast,
                        quantities = TRUE,
                        fills = list(alpha = 0.6, col = c("#4C78A8", "#F58518")),
                        edges = TRUE,
                        labels = list(font = 2),
                        main = "Bcl6 + Creb overlap (Fasted)"
                      )

creb_all_venn
creb_fed_venn
creb_fast_venn

# Other ChIP analysis - distribution, TSSdistance, annoTypes-----------------------------------------------------

# fed_peak = GenomicRanges::GRangesList(bcl6_fed, creb_peaks_fed)
# fast_peak = GenomicRanges::GRangesList(bcl6_fast, creb_peaks_fast) 
# names(fed_peak) <- c("Bcl6", "Creb")
# names(fast_peak) <- c("Bcl6", "Creb")

bcl6_fed_cov <- covplot(bcl6_fed) + labs(title = "ChIP Peaks - Bcl6 - fed")
bcl6_fast_cov <- covplot(bcl6_fast) + labs(title = "ChIP Peaks - Bcl6 - fasted")
creb_fed_cov <- covplot(creb_peaks_fed) + labs(title = "ChIP Peaks - Creb - fed")
creb_fast_cov <- covplot(creb_peaks_fast) + labs(title = "ChIP Peaks - Creb - fasted")

bcl6_fast_cov
bcl6_fed_cov
creb_fast_cov
creb_fed_cov

# this was the reference for the below figures
# plotDistToTSS(peak_anno_fast)

df1 <- as.data.frame(peak_anno_fast)
df2 <- as.data.frame(peak_anno_fed)
df3 <- as.data.frame(creb_peak_anno_fast)
df4 <- as.data.frame(creb_peak_anno_fed)


h1 <- ggplot(df1, aes(x = distanceToTSS)) +
        geom_histogram(bins = 1000) +
        coord_cartesian(xlim = c(-100000, 100000)) +
        labs(title = "Distance to TSS in Bcl6 - fasted") +
        theme_minimal()
  

h2 <- ggplot(df2, aes(x = distanceToTSS)) +
        geom_histogram(bins = 1000) +
        coord_cartesian(xlim = c(-100000, 100000)) +
        labs(title = "Distance to TSS in Bcl6 - fed") +
        theme_minimal()

h3 <- ggplot(df3, aes(x = distanceToTSS)) +
        geom_histogram(bins = 1000) +
        coord_cartesian(xlim = c(-100000, 100000)) +
        labs(title = "Distance to TSS in Creb - fast") +
        theme_minimal()

h4 <- ggplot(df4, aes(x = distanceToTSS)) +
        geom_histogram(bins = 1000) +
        coord_cartesian(xlim = c(-100000, 100000)) +
        labs(title = "Distance to TSS in Creb - fed") +
        theme_minimal()

TSSplots <- (h1 + h2) / (h3 + h4)
TSSplots

p1 <- plotAnnoBar(peak_anno_fast) + labs(title = "Feature Distribution - Bcl6 - fasted")
p2 <- plotAnnoBar(peak_anno_fed) + labs(title = "Feature Distribution - Bcl6 - fed")
p3 <- plotAnnoBar(creb_peak_anno_fast) + labs(title = "Feature Distribution - Creb - fasted")
p4 <- plotAnnoBar(creb_peak_anno_fed) + labs(title = "Feature Distribution - Creb - fed")

Annoplots <- (p1 + p2) / (p3 + p4)  + plot_layout(guides = "collect") & theme(legend.position = "right")
Annoplots

# Creb GO -----------------------------------------------------------------

ego_creb_fed1 <- enrichGO(
  gene = bcl6_creb_overlap_fed,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
ego_creb_fed2 <- enrichGO(
  gene = bcl6_creb_overlap_fed,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "MF", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
ego_creb_fed3 <- enrichGO(
  gene = bcl6_creb_overlap_fed,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "CC", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

crebfedgo1 <- dotplot(ego_creb_fed1, title = "Bcl6 + Creb - fed (Biological Process)") + scale_y_discrete(labels = function(x) str_wrap(x, width = 60))
crebfedgo2 <- dotplot(ego_creb_fed2, title = "Bcl6 + Creb - fed (Molecular Function)") + scale_y_discrete(labels = function(x) str_wrap(x, width = 60))
crebfedgo3 <- dotplot(ego_creb_fed3, title = "Bcl6 + Creb - fed (Cellular Component)") + scale_y_discrete(labels = function(x) str_wrap(x, width = 60))
crebfedgoplot <- crebfedgo1 / crebfedgo2 / crebfedgo3 + plot_layout(guides = "collect") & theme(legend.position = "right", axis.text.y = element_text(size = 8))
crebfedgoplot

ego_creb_fast1 <- enrichGO(
  gene = bcl6_creb_overlap_fast,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
ego_creb_fast2 <- enrichGO(
  gene = bcl6_creb_overlap_fast,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "MF", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
ego_creb_fast3 <- enrichGO(
  gene = bcl6_creb_overlap_fast,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "CC", #"BP" is most common
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

crebfastgo1 <- dotplot(ego_creb_fast1, title = "Bcl6 + Creb - fasted (Biological Process)") + scale_y_discrete(labels = function(x) str_wrap(x, width = 60))
crebfastgo2 <- dotplot(ego_creb_fast2, title = "Bcl6 + Creb - fasted (Molecular Function)") + scale_y_discrete(labels = function(x) str_wrap(x, width = 60))
crebfastgo3 <- dotplot(ego_creb_fast3, title = "Bcl6 + Creb - fasted (Cellular Component)") + scale_y_discrete(labels = function(x) str_wrap(x, width = 60))
crebfastgoplot <- crebfastgo1 / crebfastgo2 / crebfastgo3 + plot_layout(guides = "collect") & theme(legend.position = "right", axis.text.y = element_text(size = 8))
crebfastgoplot



# Bcl6 RNASeq -------------------------------------------------------------

# warnings = at the end, many genes have 0 expression and NA logfold change
bcl6_seq_data <- read_excel("GSE118787_Differential-Expression-Male-Livers.xlsx", 
                            sheet = "Fed WTvLKO") %>% 
                 clean_names() %>%
                 rename("symbol" = "name",
                        "Log_FC" = "wt_vs_lko_log2_fold_change",
                        "pval" = "wt_vs_lko_p_value",
                        "pval_adj" = "wt_vs_lko_adj_p_value"
                       ) %>%
                 mutate(ur_dr = case_when(
                   is.na(Log_FC) ~ NA_character_,
                   Log_FC > 0 ~ "up",
                   Log_FC < 0 ~ "down",
                   Log_FC == 0 ~ "neutral",
                   TRUE ~ NA_character_
                 ))
                        



# Filtering Bcl6 RNASeq results -------------------------------------------

# create new df with minimum needed values and only include genes with same directionality in bcl6 and gcgr ko
KO_filtered_genelist <- bcl6_seq_data %>%
  select(symbol, ur_dr, Log_FC, pval, pval_adj) %>%
  inner_join(
    gcgr_KO_data %>%
      select(symbol, ur_dr, Log_FC, pval, pval_adj) %>%
      rename(
        KO_Log_FC = Log_FC,
        KO_pval = pval,
        KO_pval_adj = pval_adj,
        KO_ur_dr = ur_dr
      ), 
    by = "symbol"
  ) %>%
  filter (ur_dr == KO_ur_dr) %>%
  select(symbol, ur_dr, Log_FC, pval, pval_adj, KO_ur_dr, KO_Log_FC, KO_pval, KO_pval_adj)

# visualize
ko_volcano <- volcano_plotter(KO_filtered_genelist, 0.05, 0.6, 260, 10, "DEGs in Bcl6 LKO with same directionality in Gcgr -/- KO", rownames = FALSE)
ko_volcano


# identical filtering/processing for gcgr mabkd
KD_filtered_genelist <- bcl6_seq_data %>%
  select(symbol, ur_dr, Log_FC, pval, pval_adj) %>%
  inner_join(
    gcgr_mAbKD_data %>%
      select(symbol, ur_dr, Log_FC, pval, pval_adj) %>%
      rename(
        KD_Log_FC = Log_FC,
        KD_pval = pval,
        KD_pval_adj = pval_adj,
        KD_ur_dr = ur_dr
      ), 
    by = "symbol"
  ) %>%
  filter (ur_dr == KD_ur_dr) %>%
  select(symbol, ur_dr, Log_FC, pval, pval_adj, KD_ur_dr, KD_Log_FC, KD_pval, KD_pval_adj)

kd_volcano <- volcano_plotter(KD_filtered_genelist, 0.05, 0.6, 270, 9, "DEGs in Bcl6 LKO with same directionality in Gcgr mAb KD", rownames = FALSE)
kd_volcano

# save plots to PDF
ggsave(
  filename = "ko_volcano_plot.pdf",
  plot = ko_volcano,
  width = 7,
  height = 5
)

ggsave(
  filename = "kd_volcano_plot.pdf",
  plot = kd_volcano,
  width = 7,
  height = 5
)

# compare gene symbols to ChIP data
KO_filtered_genelist <- KO_filtered_genelist %>%
  mutate(bcl6_chip_assoc = if_else(symbol %in% bcl6_symbols_fed, "Y", "N")) %>%
  mutate(creb_chip_assoc = if_else(symbol %in% creb_symbols_fed, "Y", "N"))

KD_filtered_genelist <- KD_filtered_genelist %>%
  mutate(bcl6_chip_assoc = if_else(symbol %in% bcl6_symbols_fed, "Y", "N")) %>%
  mutate(creb_chip_assoc = if_else(symbol %in% creb_symbols_fed, "Y", "N"))

# look at list of genes with peaks associated with creb and bcl6
shortlist_KO <- KO_filtered_genelist %>%
                  filter(bcl6_chip_assoc == "Y",
                  creb_chip_assoc == "Y")

shortlist_KD <- KD_filtered_genelist %>%
                  filter(bcl6_chip_assoc == "Y",
                  creb_chip_assoc == "Y")

# export complete filtered lists and shortlists to excel
export_list <- list(
  "Bcl6 DEGs also in Gcgr KO" = KO_filtered_genelist,
  "Bcl6 DEGs also in Gcgr mAb KD" = KD_filtered_genelist,
  "Bcl6Gcgr KO w assoc ChiP peaks" = shortlist_KO,
  "Bcl6Gcgr KD w assoc ChIP peaks" = shortlist_KD
)

write.xlsx(export_list, file = "Bcl6_Gcgr_RNASeq_Overlap.xlsx")
    
# Recreational Metabolome Prediction --------------------------------------

# chip peaks for this analysis
bcl6_symbols_fed

# need to extract DEGs from rnaseq
bcl6_degs <- bcl6_seq_data %>%
  filter(pval_adj < 0.05, abs(Log_FC) > 1)

bcl6_deg_symbols <- bcl6_degs$symbol

# intersect with chip targets
direct_targets <- intersect(bcl6_deg_symbols, bcl6_symbols_fed)

# consistent naming of direct vs all targets
deg_all <- bcl6_degs$symbol
deg_direct <- direct_targets
deg_indirect <- setdiff(deg_all, deg_direct)
bcl6_degs_direct <- bcl6_degs %>% filter(symbol %in% deg_direct)
bcl6_degs_indirect <- bcl6_degs %>% filter(symbol %in% deg_indirect)

# pull effect size and gene symbol
gene_rank_all <- bcl6_degs$Log_FC
names(gene_rank_all) <- bcl6_degs$symbol
gene_rank_all <- sort(gene_rank_all, decreasing = TRUE)
gene_rank_direct <- bcl6_degs_direct$Log_FC
names(gene_rank_direct) <- bcl6_degs_direct$symbol
gene_rank_direct <- sort(gene_rank_direct, decreasing = TRUE)



# Plot Calls + Merges --------------------------------------------------------------
kdkovenn
vp1
vp2
kogoplot
kdgoplot
creb_all_venn
creb_fed_venn
creb_fast_venn
bcl6_fast_cov
bcl6_fed_cov
creb_fast_cov
creb_fed_cov
TSSplots
Annoplots
crebfedgoplot
crebfastgoplot
ko_volcano
kd_volcano


