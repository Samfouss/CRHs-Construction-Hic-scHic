# Chargement des librairies
source("./scripts/load_save_data.R")
source("./scripts/crhs_creation_with_cells/create_graph_from_cells.R")

# Présentation des clusters
table(cluster_with_imp$cluster)

nb_cells = 250
cell_dim = 562
# Conversion des valeurs de la matrice de densité en 0 et 1.
# 1 si la valeur est sipérieure à 0 et 0 sinon
MCMCImpute_result_Impute_SZ_ <- MCMCImpute_result$Impute_SZ > 0

# Boucle sur les données imputées pour construire les matrices d'incidence
clusters = unique(cluster_with_imp$cluster)
nb_clu = length(clusters)

clu_chrs = list()
for (clu in clusters) {
  
  # On recupère ici les cellules qui font parties du cluster en cours
  class = cluster_with_imp$cluster[cluster_with_imp$cluster==clu]
  # On recupère ici les numero des cellules
  elems = as.numeric(str_sub(names(class), 5, 7))
  
  # On boucle sur les cellules du cluster courant
  for (cell in elems) {
    
    # Ce seuil nous permet de detecter les concetrations que l'on considère comme du bruit ou comme un vrai contact dans les données retournées par HicImput. Ici on retient par la moyenne des concentrations dans une cellule données
    seuil = mean(MCMCImpute_result$Impute_SZ[, cell][MCMCImpute_result$Impute_SZ[, cell]!=0])
    temp <- as.numeric(MCMCImpute_result$Impute_SZ[, cell] > seuil)
    
    mat_row_name = paste0("Loci", seq_len(cell_dim))
    incidence_mat <- matrix(
      0, 
      nrow = cell_dim, 
      ncol = cell_dim,
      dimnames = list(mat_row_name, mat_row_name)
    )

    incidence_mat[upper.tri(incidence_mat, diag = FALSE)] <- temp
    incidence_mat <- incidence_mat + t(incidence_mat)
    # On verifie bien que la matrice est symétrique
    # isSymmetric(incidence_mat)
    
    # Initialisation du premier cluster
    print("Cluster")
      
    all_net_result <- create_bip_clust_graph_from_cell(
      incidence_mat, 
      scHic_promoters_ids, 
      cell
    )
  }
  
  clu_chrs[[length(clu_chrs)+1]] <- all_net_result
}

names(clu_chrs) <- str_c("cluster", clusters)
clu_chrs

save(clu_chrs, file = "rdata/clu_chrs_result.rda")

