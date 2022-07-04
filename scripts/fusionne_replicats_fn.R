
fusion_comp <- function(structure_net_comp1, structure_net_comp2, edge_overlap_result){
  
  #' @description Fusionne deux CRHs appariés dans des réplicats différents
  #' @param structure_net_comp1 Composantes du réplicat 1 retournées par la fonction components
  #' @param structure_net_comp2 Composantes du réplicat 2 retournées par la fonction components
  #' @param edge_overlap_result Vecteur de longueur 2 avec l'indice de la composante du réplicat 1 et de celle du réplicat 2 à fusionner
  
  results = list()
  n_block <- length(structure_net_comp1)
  
  for (b in seq_len(n_block)) {
    paires.crh <- edge_overlap_result[[b]]$chev_edge_comp
    for (p in seq_len(ncol(paires.crh))) {
      
      # Matrice d'adjacence du réplicat 1
      dist_bin1 <- structure_net_comp1[[b]]$dist_bin
      # Matrice d'adjacence du réplicat 2
      dist_bin2 <- structure_net_comp2[[b]]$dist_bin
      
      ens1 = names(structure_net_comp1[[b]]$membership_bip)[structure_net_comp1[[b]]$membership_bip==paires.crh[1, p]]
      ens2 = names(structure_net_comp2[[b]]$membership_bip)[structure_net_comp2[[b]]$membership_bip==paires.crh[2, p]]
      noms = union(ens1,ens2)
      adjcomp = matrix(0,length(noms),length(noms))
      dimnames(adjcomp) = list(noms,noms)
      
      adjcomp[ens1,ens1] = dist_bin1[ens1,ens1]
      adjcomp[ens2,ens2] = adjcomp[ens2,ens2]|dist_bin2[ens2,ens2]
      
      results[[length(results)+1]] <- adjcomp
      
    }
  }
  
  names(results) <- str_c("adjcomp", seq_len(n_block))
  results
  
}


