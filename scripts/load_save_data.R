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
library(data.table)
#knitr::knit_hooks$set(webgl = hook_webgl)
# On Charge les IDs des promoters deja sauvés
load("rdata/promoters_ids.rda")
# On charge ici dans l'environnement toutes les données des structures.
# structure_data contient 250 éléments. Chaque élément contient une paire de structure (paire matchée dans le dossier data/pairs sur gitHub)
load("rdata/structure_data.rda")

structure_1 <- structure_data$paired_1$structure_1%>%select(-c(ends_with("_c"), "replica"))
structure_2 <- structure_data$paired_1$structure_2%>%select(-c(ends_with("_c"), "replica"))

structure_3 <- structure_data$paired_2$structure_1%>%select(-c(ends_with("_c"), "replica"))
structure_4 <- structure_data$paired_2$structure_2%>%select(-c(ends_with("_c"), "replica"))

# Création d'un identifiant pour chaque lignes
structure_1$ID <- paste0("B", sprintf("%02d", structure_1$X4), sprintf("%04d", 1:nrow(structure_1)))
structure_2$ID <- paste0("B", sprintf("%02d", structure_2$X4), sprintf("%04d", 1:nrow(structure_2)))


