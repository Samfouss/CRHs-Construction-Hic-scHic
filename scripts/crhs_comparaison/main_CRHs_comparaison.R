# Chargement des librairies
library(tidyverse)
library(igraph)

# Etape 1 : Chargement des CRHs obtenus a partir des structures
load("rdata/all_rda_data/all_net_result_complex_.rda")
#load("rdata/all_net_result_complex.rda")

# Les résultats de quel cluster ?
nb_clusters = 6

# Etape 2 : Chargement des CRHs obtenus à partir des clusters
load(paste0("rdata/all_rda_data/clu_chrs_result_", nb_clusters, ".rda"))

# Etape 3 : fonction permettant de faire la réduction de matrices
source("scripts/crhs_comparaison/degenerationMatrix_fn.R")

# Etape 4 : fonction permettant de faire la comparaison de CRHs
source("scripts/crhs_comparaison/compute_comparaison.R")

##################  Comparaison des CRHs ####################

crhs_comparation_res = compute_comparaison(all_net_result_complex_, clu_chrs_result, make_degeneration = FALSE)

########### Sauvegarde des données ###########
save(crhs_comparation_res, file = paste0("rdata/all_rda_data/crhs_comparation_res_", nb_clusters, ".rda"))

