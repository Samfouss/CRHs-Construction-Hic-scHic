#################### Etapes de Création de CRHs à partir des structures et des clusters de cellules ####################

# Chargement des librairies

library(stats)
library(stringr)
library(FactoMineR)
library(factoextra)
library(fpc)
library(cluster)
library(stats)
library(ggpubr)
library(igraph)
library(tidyverse)

### Création des CRHs à partir des structures

source("scripts/crhs_creation_with_bins/main_crh_creation_with_structure.R")

### Création des CRHs à partir des clusters de cellules

# Les résultats de quel cluster ?
nb_clusters = 7

source("scripts/crhs_creation_with_cells/main_CRHs_creation_with_cells.R")

### Comparaison des CRHs 

source("scripts/crhs_comparaison/main_CRHs_comparaison.R")

