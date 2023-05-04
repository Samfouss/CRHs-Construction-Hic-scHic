

edge_identity_overlap_from_cell <- function(scHic_net_comp1, scHic_net_comp2){

  # On recupere ici le nombre de CRHs des deux compposantes
  nb_crh1 = length(scHic_net_comp1$crhs)
  nb_crh2 = length(scHic_net_comp2$crhs)
  
  chevauche_1_2 = matrix(NA, nb_crh1, nb_crh2)
  
  if(sum(scHic_net_comp1$crhs[[1]]$mat_incidence)==-1){
    
    results <- list(
      "chevauche_1_2_perc" = 0,
      "chev_edge_comp1" = 0,
      "chev_edge_comp2" = 0,
      "chev_edge_comp" = 0
    )
    
  }else if(sum(scHic_net_comp2$crhs[[1]]$mat_incidence)==-1){
    
    results <- list(
      "chevauche_1_2_perc" = 0,
      "chev_edge_comp1" = 0,
      "chev_edge_comp2" = 0,
      "chev_edge_comp" = 0
    )
    
  }else{
    
    for (i in seq_len(nb_crh1)) {
      
      # On recupère ici la matrice d'incidence
      mat1 = scHic_net_comp1$crhs[[i]]$mat_incidence
      
      for (j in seq_len(nb_crh2)) {
        
        # On recupère ici la matrice d'incidence
        mat2 = scHic_net_comp2$crhs[[j]]$mat_incidence
        
        # Ici je recupère les noms des lignes et colonnes
        rowmat = union(rownames(mat1), rownames(mat2))
        colmat = union(colnames(mat1), colnames(mat2))
        
        # Je crée ici des matrices qui doivent recevoir les données des matrices redimenssionnées. J'en crée 2, une pour mat1 et une seconde pour mat2
        mat_merge_1 <- matrix(
          0,
          ncol=length(colmat),
          nrow=length(rowmat),
          dimnames=list(rowmat, colmat)
        )
        mat_merge_2 <- mat_merge_1
        
        # On récupère ici les indexes des cellules contenant les données de depart afin de n'ajouter que ces derniers dans les matrices redimensionnées. Le reste est rempli par des zéros.
        indxA <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat1), colnames(mat1), FUN=paste)
        indxB <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat2), colnames(mat2), FUN=paste)
        # On récupère ici les données dans les matrices redimenssionnées. Avec les noms des lignes et colonnes (indexes), on est sur de recuprerer les bonnes données
        mat_merge_1[indxA] <- mat1
        mat_merge_2[indxB] <- mat2
        
        # On calcul les taux d'identité
        chevauche_1_2[i,j] = sum(mat_merge_1&mat_merge_2, na.rm = TRUE)/sum(mat_merge_1|mat_merge_2, na.rm = TRUE)
        
      }
      
    }
    
    c1.list = apply(chevauche_1_2,1,which.max)
    c2.list = apply(chevauche_1_2,2,which.max)
    
    chev_edge_comp1 = rbind(1:nb_crh1, c1.list)
    # Composantes du réplicat 2 et leur composante chevauchant le plus dans le réplicat 1
    chev_edge_comp2 = rbind(1:nb_crh2, c2.list)
    
    chev_edge_comp = matrix(NA, 2, ncol = min(ncol(chev_edge_comp1), ncol(chev_edge_comp2)))
    
    for (m1 in seq_len(ncol(chev_edge_comp1))) {
      for (m2 in seq_len(ncol(chev_edge_comp2))) {
        if(all(chev_edge_comp1[, m1]==rev(chev_edge_comp2[, m2]))){
          chev_edge_comp[, sum(!is.na(chev_edge_comp[1, ])) + 1] <- chev_edge_comp1[, m1]
        } 
      }
    }
    
    chev_edge_comp <- matrix(chev_edge_comp[!is.na(chev_edge_comp)], 2, ncol = length(chev_edge_comp[!is.na(chev_edge_comp)])/2)
    
    
    results <- list(
      "chevauche_1_2_perc" = chevauche_1_2,
      "chev_edge_comp1" = chev_edge_comp1,
      "chev_edge_comp2" = chev_edge_comp2,
      "chev_edge_comp" = chev_edge_comp
    )
    
  }
  
  results
  
}
