

#' Cette fonction permet d'appréhender le taux de chevauchement entre les CRHs
#'
#' @description
#' #' `edge_identity_overlap` retourne la similarité des CRHs en terme de nombre d'arrêts. L'idée est de partir de la matrice binaire des distance entre les différentes billes et de compter exactement le nombre de billes qui concident en terme de nombre d'arrêtes (correspond a 1 dans la matrice binaire) sur le nombre de billes total

edge_identity_overlap <- function(structure_net_comp1, structure_net_comp2, seuil=0){
  
  # On recupère les blocs à parcourir
  n_block = max(
    length(structure_net_comp1), 
    length(structure_net_comp2)
  )
  # Initialisation de la liste à recevoir les résultats
  results = list()
  
  for (b in seq_len(n_block)){
    
    # On recupere ici le nombre de CRHs des deux compposantes
    nb_crh1 = length(structure_net_comp1[[b]]$crhs)
    nb_crh2 = length(structure_net_comp2[[b]]$crhs)
    
    chevauche_1_2 = matrix(0, nb_crh1, nb_crh2)
    
    if(sum(structure_net_comp1[[b]]$crhs[[1]]$mat_incidence)==-1 || sum(structure_net_comp2[[b]]$crhs[[1]]$mat_incidence)==-1){
      
      results[[length(results)+1]] <- list(
        "chevauche_1_2_perc" = 0,
        "chev_edge_comp1" = 0,
        "chev_edge_comp2" = 0,
        "chev_edge_comp" = 0
      )
      
    }else{
      
      for (i in seq_len(nb_crh1)) {
        
        # On recupère ici la matrice d'incidence
        mat1 = structure_net_comp1[[b]]$crhs[[i]]$mat_incidence
        
        for (j in seq_len(nb_crh2)) {
          
          # On recupère ici la matrice d'incidence
          mat2 = structure_net_comp2[[b]]$crhs[[j]]$mat_incidence
          
          
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
      chevauche = chevauche_1_2
      chevauche_1_2[chevauche_1_2<=seuil]=NA
      
      if(sum(chevauche_1_2, na.rm = TRUE) == 0){
        
        results[[length(results)+1]] <- list(
          "chevauche_1_2_perc" = 0,
          "chev_edge_comp1" = 0,
          "chev_edge_comp2" = 0,
          "chev_edge_comp" = 0
        )
        
      }else{
        
        # chevauche_1_2[rowSums(chevauche_1_2, na.rm = TRUE)>0, , drop=FALSE]
        c1.list = unlist(apply(chevauche_1_2,1,which.max))
        c2.list = unlist(apply(chevauche_1_2,2,which.max))
        
        # Lorsque la colonne vaut -1, cela voudrait dire qu'il y a pas de chevauchement avec le seuil choisi sur la ligne
        chev_edge_comp1 = rbind(1:nb_crh1, -1)
        chev_edge_comp1[2, rowSums(chevauche_1_2, na.rm = TRUE)>0] <- c1.list
        # Je retir les colonnes où la ligne 2 vaut -1
        chev_edge_comp1 = chev_edge_comp1[, which(chev_edge_comp1[2, ]!=-1), drop=FALSE]
        # Composantes du réplicat 2 et leur composante chevauchant le plus dans le réplicat 1
        #       # Lorsque la colonne vaut -1, cela voudrait dire qu'il y a pas de chevauchement avec le seuil choisi sur la colonne
        chev_edge_comp2 = rbind(1:nb_crh2, -1)
        chev_edge_comp2[2, colSums(chevauche_1_2, na.rm = TRUE)>0] <- c2.list
        # Je retir les colonnes où la ligne 2 vaut -1
        chev_edge_comp2 = chev_edge_comp2[, which(chev_edge_comp2[2, ]!=-1), drop=FALSE]
        
        chev_edge_comp = matrix(NA, 2, ncol = min(ncol(chev_edge_comp1), ncol(chev_edge_comp2)))
      # print(ncol(chev_edge_comp1))
      # print(ncol(chev_edge_comp2))
        for (m1 in seq_len(ncol(chev_edge_comp1))) {
          for (m2 in seq_len(ncol(chev_edge_comp2))) {
            if(all(chev_edge_comp1[, m1]==rev(chev_edge_comp2[, m2]))){
              chev_edge_comp[, sum(!is.na(chev_edge_comp[1, ])) + 1] <- chev_edge_comp1[, m1]
            } 
          }
        }
      
        chev_edge_comp <- matrix(chev_edge_comp[!is.na(chev_edge_comp)], 2, ncol = length(chev_edge_comp[!is.na(chev_edge_comp)])/2)
        
        results[[length(results)+1]] <- list(
          "chevauche_1_2_perc" = chevauche,
          "chev_edge_comp1" = chev_edge_comp1,
          "chev_edge_comp2" = chev_edge_comp2,
          "chev_edge_comp" = chev_edge_comp
        )
      }
      
    }
    
  }
  
  names(results) <- str_c("edge_overlap", seq_len(length(results)))
  results
  
}


# res = edge_identity_overlap2(net1, net2)



