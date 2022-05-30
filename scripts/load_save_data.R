# Chargement des librairies
library(rmarkdown)
library(rgl)
library(magick)
library(randomcoloR)
library(readr)
library(tidyverse)
library(igraph)
library(ggpubr)
library(RColorBrewer)
library(finalfit)
library(knitr)
library(plotly)
knitr::knit_hooks$set(webgl = hook_webgl)

# Initialisation des paramètres
rep_number = 4

# Importation des données
for (j in 1:rep_number) {
  data_url = paste0("https://raw.githubusercontent.com/fmusella/In-silico_Hi-C_GAM_SPRITE/main/data/structures/structure_", j, ".txt")
  data_name = paste("structure_", j, sep="")
  
  assign(
    data_name, 
    read_table(
      data_url, 
      col_names = FALSE
    )
  )
}

# Création d'un identifiant pour chaque lignes
structure_1$ID <- paste0("B", structure_1$X4, sprintf("%03d", 1:nrow(structure_1)))
structure_2$ID <- paste0("B", structure_2$X4, sprintf("%03d", 1:nrow(structure_2)))
structure_3$ID <- paste0("B", structure_3$X4, sprintf("%03d", 1:nrow(structure_3)))
structure_4$ID <- paste0("B", structure_4$X4, sprintf("%03d", 1:nrow(structure_4)))


