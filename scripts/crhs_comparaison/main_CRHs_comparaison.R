# Ce fichier permet de faire la comparaison et stocke les résultats dans le fichier crhs_comparation_res

# Etape 1 : Chargement des CRHs obtenus a partir des structures dites complexes et après elimination des réseaux redondants
load("rdata/all_rda_data/all_net_result_complex2_.rda")
#load("rdata/all_rda_data/all_net_result_complex.rda")

# Etape 2 : Chargement des CRHs obtenus à partir des clusters
load(paste0("rdata/all_rda_data/clu_chrs_result_", nb_clusters, ".rda"))

# Etape 3 : fonction permettant de faire la réduction de matrices
source("scripts/crhs_comparaison/degenerationMatrix_fn.R")

# Etape 4 : fonction permettant de faire la comparaison de CRHs
source("scripts/crhs_comparaison/compute_comparaison.R")

##################  Comparaison des CRHs ####################

crhs_comparation_res = compute_comparaison(all_net_result_complex2_, clu_chrs_result, make_degeneration = FALSE)

########### Sauvegarde des données ###########
save(crhs_comparation_res, file = paste0("rdata/all_rda_data/crhs_comparation_res_", nb_clusters, ".rda"))

