
# Etape 1 : Chargement des CRHs obtenus a partir des structures
load("rdata/all_net_result.rda")

# Etape 2 : Chargement des CRHs obtenus à partir des clusters
load("rdata/clu_chrs_result.rda")

# Etape 3 : fonction permettant de faire la réduction de matrices
source("scripts/crhs_comparaison/degenerationMatrix_fn.R")

# Etape 4 : fonction permettant de faire la comparaison de CRHs
source("scripts/crhs_comparaison/compute_comparaison.R")

##################  Comparaison des CRHs ####################

crhs_comparation_res = compute_comparaison(all_net_result, clu_chrs_result)




