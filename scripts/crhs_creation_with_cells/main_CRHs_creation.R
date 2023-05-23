# Chargement des librairies
library(igraph)
library(dplyr)
library(stringr)
library(stats)
library(data.table)
library(tidyverse)

# Etape 1 : Transformation des données des matrices scHic
# source("scripts/crhs_creation_with_cells/getMatrixUperDataAndSave.R")

# Etape 2 :  Création de clusters
# source("scripts/crhs_creation_with_cells/cells_clustering.R")

# Etape 3 : Transformation des sorties de Juicer en matrice d'incidence pour la création de graphs
# source("scripts/crhs_creation_with_cells/fromJuicerOutputFileToMatrix.R")

# Etape 4 : Chargement de la fonction pour la création des graphs
source("scripts/crhs_creation_with_cells/create_graph_from_cells.R")

# Chargement des promoters
load("rdata/scHic_promoters_ids.rda")

# Chargement des données sur les clusters
load("rdata/merge_loops_clus.rda")

clu_chrs_result = list()
for (clus in seq_len(length(merge_loops_clus))) {
 
  net = create_bip_clust_graph_from_cell(merge_loops_clus[[clus]], scHic_promoters_ids, clus)
  clu_chrs_result[[length(clu_chrs_result)+1]] <- net
}

names(clu_chrs_result) <- str_c("cluster", seq_len(length(merge_loops_clus)))

############## Distribution des crhs dans les 5 clusters ##########
ncrhs <- matrix(
  0,
  nrow = length(clu_chrs_result),
  ncol = 1,
)
for (res in seq_len(length(clu_chrs_result))) {
  # sum(fireCaller_result[[res]]$FIRE_output$Mousse_cells_clus_12_indicator)
  ncrhs[res, 1] <- length(clu_chrs_result[[res]]$crhs)
}

summary(ncrhs[, 1])

########### Sauvegarde des données ###########
save(clu_chrs_result, file = "rdata/clu_chrs_result.rda")

