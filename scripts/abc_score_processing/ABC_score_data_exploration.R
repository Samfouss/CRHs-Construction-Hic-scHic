library(tidyverse)
library(data.table)
library(vroom)

endo_cell_data <- read_tsv("rdata/SCORE_ABC_RESULTS/ENDO_CELL/EnhancerPredictions.txt", show_col_types = FALSE)

# To charge large dataset
endo_cell_data_ <- fread("rdata/SCORE_ABC_RESULTS/ENDO_CELL/EnhancerPredictionsAllPutative.txt")

# endo_cell_data_ <- vroom("rdata/SCORE_ABC_RESULTS/ENDO_CELL/EnhancerPredictionsAllPutative.txt", delim = "\t")

seuils <- seq(0, 0.10, 0.001)
contact_G_H_dist <- matrix(
  NA, 
  nrow = length(seuils), 
  ncol = 9,
  dimnames = list(
    c(),
    c("threshold", "Reject", "Keep", "min_TarGene", "frst_Quar_TarGene", "median_TarGene", "mean_TarGene", "third_Quar_TarGene", "max_TarGene"))
)

for (i in seq_len(length(seuils))) {
  contact_G_H_dist[i, 1] <- seuils[i]
  contact_G_H_dist[i, 2] <- table(endo_cell_data$ABC.Score>seuils[i])[1]
  contact_G_H_dist[i, 3] <- table(endo_cell_data$ABC.Score>seuils[i])[2]
  
  summarise <- data.frame(
    table(endo_cell_data%>%filter(ABC.Score>seuils[i])%>%select(TargetGene))
  )%>%
    summarise(
      min_TarGene = min(Freq, na.rm = TRUE),
      frst_Quar_TarGene = quantile(Freq, 0.25, na.rm = FALSE),
      median_TarGene = quantile(Freq, 0.5, na.rm = FALSE),
      mean_TarGene = mean(Freq, na.rm = TRUE),
      third_Quar_TarGene = quantile(Freq, 0.75, na.rm = FALSE),
      max_TarGene = max(Freq, na.rm = TRUE)
    )
  
  contact_G_H_dist[i, 4] <- summarise$min_TarGene
  contact_G_H_dist[i, 5] <- summarise$frst_Quar_TarGene
  contact_G_H_dist[i, 6] <- summarise$median_TarGene
  contact_G_H_dist[i, 7] <- summarise$mean_TarGene
  contact_G_H_dist[i, 8] <- summarise$third_Quar_TarGene
  contact_G_H_dist[i, 9] <- summarise$max_TarGene
}

hist(contact_G_H_dist[, 2])

# Visualiser les données sur les genes retenus et rejetés
contact_G_H_dist <- data.frame(contact_G_H_dist)
ggplot(contact_G_H_dist) + geom_col(aes(x = threshold, y = Reject))
ggplot(contact_G_H_dist) + geom_col(aes(x = threshold, y = Keep))
ggplot(contact_G_H_dist) + geom_col(aes(x = threshold, y = mean_TarGene))



