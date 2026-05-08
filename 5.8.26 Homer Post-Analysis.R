
# 0. Initialize -----------------------------------------------------------
pacman::p_load(
  ggplot2,
  dplyr,
  reshape2,
  try.bioconductor = TRUE
)

# 1. Read File ------------------------------------------------------------

union <- read_tsv("peak_tag_matrix_union.txt")
intersection <- read_tsv("peak_tag_matrix_intersection.txt")
KO <- read_tsv("peak_tag_matrix_KO.txt")
KD <- read_tsv("peak_tag_matrix_KD.txt")

new_names <- c("tag_LKO_GCG_Input", "tag_LKO_Veh_Input", "tag_P8249_FLKO_GCG", "tag_P8250_FLKO_GCG", "tag_P8251_FWT_GCG", "tag_P8266_FWT_GCG", "tag_P8268_FLKO_Veh", "tag_P8289_FLKO_Veh", "tag_P8292_FWT_Veh", "tag_P8297_FWT_Veh", "tag_WT_GCG_Input", "tag_WT_Veh_Input")

names(union)[20:31] <- new_names
names(intersection)[20:31] <- new_names
names(KO)[20:31] <- new_names
names(KD)[20:31] <- new_names


# 2. Make boxplot ---------------------------------------------------------

start <- 22
end <- 29

union_data <- union[, start:end]
intersection_data <- union[, start:end]
KO_data <- union[, start:end]
KD_data <- union[, start:end]


union_long <- melt(union_data)
intersection_long <- melt(intersection_data)
KO_long <- melt(KO_data)
KD_long <- melt(KD_data)

union_plot <- ggplot(union_long, aes(x = variable, y = value)) +
  geom_boxplot(
    # outlier.shape = NA
    ) +
  # coord_cartesian(ylim = c(0, 50)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Tag counts",
       title = "union")

intersection_plot <- ggplot(intersection_long, aes(x = variable, y = value)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Tag counts",
       title = "intersection")

KO_plot <- ggplot(KO_long, aes(x = variable, y = value)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Tag counts",
       title = "KO")

KD_plot <- ggplot(KD_long, aes(x = variable, y = value)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Tag counts",
       title = "KD")

union_plot
intersection_plot
KO_plot
KD_plot
# repeat homer workflow with the new file of total pCREB and then repeat boxplots with that too