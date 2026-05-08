# Initialize --------------------------------------------------------------
pacman::p_load(
  dplyr,
  readr,
  openxslx,
  eulerr,
  try.bioconductor = TRUE
)


# 1. Import Ellen's File ----------------------------------------------------

ellen_genelist <- read_tsv("Ellen pCREB files/pCREB.fed.BCL6.mergepeaks.anno.txt")

pcreb_total <- read_tsv("Ellen pCREB files/pCREB.all.mergepeaks.txt")

# 2. Make a list of all of my identified genes ----------------------------

KD_list <- shortlist_KD$symbol
KO_list <- shortlist_KO$symbol
combined_list <- c(KO_list, KD_list)
intersection_list <- intersect(KD_list, KO_list)

# 3. Only keep rows in Ellen's genelist that correspond to my li --------

filtered_ellen_genelist <- ellen_genelist %>%
  filter(`Gene Name` %in% combined_list)

KO_filtered_ellen_genelist <- ellen_genelist %>%
  filter(`Gene Name` %in% KO_list)

KD_filtered_ellen_genelist <- ellen_genelist %>%
  filter(`Gene Name` %in% KD_list)

intersection_filtered_ellen_genelist <- ellen_genelist %>%
  filter(`Gene Name` %in% intersection_list)

length(unique(combined_list))
length(unique(KO_list))
length(unique(KD_list))
length(unique(intersection_list))


length(unique(filtered_ellen_genelist$`Gene Name`))
length(unique(KO_filtered_ellen_genelist$`Gene Name`))
length(unique(KD_filtered_ellen_genelist$`Gene Name`))
length(unique(intersection_filtered_ellen_genelist$`Gene Name`))


# 4. Export ---------------------------------------------------------------
wb <- createWorkbook()

addWorksheet(wb, "Combined Genes")
addWorksheet(wb, "Genes from Gcgr KO")
addWorksheet(wb, "Genes from Gcgr KD")
addWorksheet(wb, "Intersection Genes")
addWorksheet(wb, "All pCREB Genes")

writeData(wb, "Combined Genes", filtered_ellen_genelist)
writeData(wb, "Genes from Gcgr KO", KO_filtered_ellen_genelist)
writeData(wb, "Genes from Gcgr KD", KD_filtered_ellen_genelist)
writeData(wb, "Intersection Genes", intersection_filtered_ellen_genelist)
writeData(wb, "All pCREB Genes", pcreb_total)

saveWorkbook(wb, "Filtered Peaks.xlsx", overwrite = TRUE)



# 5. Venn diagram ---------------------------------------------------------

fit <- euler(list(
  #Ved_combined = unique(combined_list),
  Ved_KO = unique(KO_list),
  Ved_KD = unique(KD_list),
  #Filtered_combined = unique(filtered_ellen_genelist$`Gene Name`),
  Filtered_KO = unique(KO_filtered_ellen_genelist$`Gene Name`),
  Filtered_KD = unique(KD_filtered_ellen_genelist$`Gene Name`)
  #ellen_list = unique(ellen_genelist$`Gene Name`)
))

plot(fit, quantities = TRUE)

