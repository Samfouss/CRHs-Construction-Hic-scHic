
source("scripts/crhs_creation_with_cells/fill_symetric_matrix_fn.R")
# Objet venant du clustering des cellules
load("rdata/all_rda_data/cluster_result.rda")

# Trasnformation des données des clusters en matrix
cluster_matrix_result = fill_symetric_matrix(cluster_result$cluster, cellUperDiagData)
save(cluster_matrix_result, file = "rdata/all_rda_data/cluster_matrix_result.rda")


