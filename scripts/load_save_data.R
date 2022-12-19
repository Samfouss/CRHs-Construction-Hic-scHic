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
#load("rdata/structure_data.rda")
load("rdata/all_single_structure.rda")
load("rdata/all_paired_structure.rda")

structure_1 <- as_tibble(
  all_paired_structure%>%
    filter(paire == str_c(1, sprintf("%03d", 1)))%>%
    select(-c(ends_with("_c"), "paire"))%>%
    mutate(
      ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n()))
    )
)

structure_2 <- as_tibble(
  all_paired_structure%>%
    filter(paire == str_c(2, sprintf("%03d", 1)))%>%
    select(-c(ends_with("_c"), "paire"))%>%
    mutate(
      ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n()))
    )
)

structure_3 <- as_tibble(
  all_paired_structure%>%
    filter(paire == str_c(1, sprintf("%03d", 2)))%>%
    select(-c(ends_with("_c"), "paire"))%>%
    mutate(
      ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n()))
    )
)

structure_4 <- as_tibble(
  all_paired_structure%>%
    filter(paire == str_c(2, sprintf("%03d", 2)))%>%
    select(-c(ends_with("_c"), "paire"))%>%
    mutate(
      ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n()))
    )
)

