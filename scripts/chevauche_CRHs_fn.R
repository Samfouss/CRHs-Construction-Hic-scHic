

#' Cette fonction permet d'appréhender le taux de chevauchement entre les CRHs
#'
#' @description
#' #' `edge_identity_overlap` retourne la similarité des CRHs en terme de nombre d'arrêts. L'idée est de partir de la matrice binaire des distance entre les différentes billes et de compter exactement le nombre de billes qui concident en terme de nombre d'arrêtes (correspond a 1 dans la matrice binaire) sur le nombre de billes total
#'
#' @param point_1 
#' @param point_2 


structure_net_comp1 <- net1
structure_net_comp2 <- net2


edge_identity_overlap2 <- function(structure_net_comp1, structure_net_comp2){
  
  
  # On recupère les blocs à parcourir
  n_block = length(structure_net_comp1)
  # Initialisation de la liste à recevoir les résultats
  results = list()
  
  for (b in seq_len(n_block)){
    
    # On recupere ici le nombre de CRHs des deux compposantes
    nb_crh1 = length(structure_net_comp1[[b]]$crhs)
    nb_crh2 = length(structure_net_comp2[[b]]$crhs)
    
    chevauche_1_2 = matrix(NA, nb_crh1, nb_crh2)
    
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
    
    
    results[[length(results)+1]] <- list(
      "chevauche_1_2_perc" = chevauche_1_2,
      "chev_edge_comp1" = chev_edge_comp1,
      "chev_edge_comp2" = chev_edge_comp2,
      "chev_edge_comp" = chev_edge_comp
    )
    
  }
  
  names(results) <- str_c("edge_overlap", seq_len(n_block))
  results
  
}


############################### Ancien code ###########################


edge_identity_overlap <- function(structure_net_comp1, structure_net_comp2){
  
  # On recupère les blocs à parcourir
  n_block = length(structure_net_comp1)
  # Initialisation de la liste à recevoir les résultats
  results = list()
  
  for (b in seq_len(n_block)){
    
    dist_bin1 = as.matrix(structure_net_comp1[[b]]$dist_bin)
    dist_bin2 = as.matrix(structure_net_comp2[[b]]$dist_bin)
    chevauche_1_2 = matrix(
      NA,
      structure_net_comp1[[b]]$no_bip,
      structure_net_comp2[[b]]$no_bip
    )
    
    
    for (i in (1:structure_net_comp1[[b]]$no_bip)[structure_net_comp1[[b]]$csize_bip>1]){
      
      mat1 <- matrix(
        0,
        nrow(dist_bin1),
        ncol(dist_bin1),
        dimnames = dimnames(dist_bin1)
      )
      
      membership1 <- structure_net_comp1[[b]]$membership_bip
      mat1[
        row.names(dist_bin1) %in% names(membership1[membership1==i]),
        names(membership1[membership1==i])
      ] = dist_bin1[
        row.names(dist_bin1) %in% names(membership1[membership1==i]),
        names(membership1[membership1==i])
      ]
      
      for (j in (1:structure_net_comp2[[b]]$no_bip)[structure_net_comp2[[b]]$csize_bip>1]){
        
        mat2 <- matrix(
          0,
          nrow(dist_bin2),
          ncol(dist_bin2),
          dimnames = dimnames(dist_bin2)
        )
        
        membership2 <- structure_net_comp2[[b]]$membership_bip
        mat2[
          row.names(dist_bin2) %in% names(membership2[membership2==j]),
          names(membership2[membership2==j])
        ] = dist_bin2[
          row.names(dist_bin2) %in% names(membership2[membership2==j]),
          names(membership2[membership2==j])
        ]
        
        ###################### Redimensionnement des matrices 1 et 2 #####
        
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
    # Composantes du réplicat 1 et leur composante chevauchant le plus dans le réplicat 2
    chev_edge_comp1 = rbind(
      (1:structure_net_comp1[[b]]$no_bip)[structure_net_comp1[[b]]$csize_bip>1], 
      unlist(c1.list[structure_net_comp1[[b]]$csize_bip>1])
    )
    # Composantes du réplicat 2 et leur composante chevauchant le plus dans le réplicat 1
    chev_edge_comp2 = rbind(
      (1:structure_net_comp2[[b]]$no_bip)[structure_net_comp2[[b]]$csize_bip>1], 
      unlist(c2.list[structure_net_comp2[[b]]$csize_bip>1])
    )
    
    chev_edge_comp = matrix(NA, 2, ncol = min(ncol(chev_edge_comp1), ncol(chev_edge_comp2)))
    
    for (m1 in seq_len(ncol(chev_edge_comp1))) {
      for (m2 in seq_len(ncol(chev_edge_comp2))) {
        if(all(chev_edge_comp1[, m1]==rev(chev_edge_comp2[, m2]))){
          chev_edge_comp[, sum(!is.na(chev_edge_comp[1, ])) + 1] <- chev_edge_comp1[, m1]
        } 
      }
    }
    
    chev_edge_comp <- matrix(chev_edge_comp[!is.na(chev_edge_comp)], 2, ncol = length(chev_edge_comp[!is.na(chev_edge_comp)])/2)
    
    results[[length(results)+1]] <- list(
      "chevauche_1_2_perc" = chevauche_1_2,
      "chev_edge_comp1" = chev_edge_comp1,
      "chev_edge_comp2" = chev_edge_comp2,
      "chev_edge_comp" = chev_edge_comp,
      "no" = c(length(unique(c1.list)), length(unique(c2.list)))
    )
    
  }
  
  names(results) <- str_c("edge_overlap", seq_len(n_block))
  results
  
}





#' Cette fonction permet d'appréhender le taux de chevauchement entre les CRHs
#'
#' @description
#' #' `dist_func` retourne la distance entre deux points
#'
#' @param point_1 
#' @param point_2 

# chr_overlap <- function(structure_net_comp1, structure_net_comp2){
#   
#   n_block = length(structure_net_comp1) - 1
#   results = list()
#   
#   for (b in seq_len(n_block)) {
#     
#     chevauche_1_2 = matrix(NA,structure_net_comp1[[b]]$no,structure_net_comp2[[b]]$no)
#     for (i in (1:structure_net_comp1[[b]]$no)[structure_net_comp1[[b]]$csize>1]){
#       for (j in (1:structure_net_comp2[[b]]$no)[structure_net_comp2[[b]]$csize>1]){
#         ens1 = names(structure_net_comp1[[b]]$membership)[structure_net_comp1[[b]]$membership==i]
#         ens2 = names(structure_net_comp2[[b]]$membership)[structure_net_comp2[[b]]$membership==j]
#         chevauche_1_2[i,j] = length(intersect(ens1,ens2)) /length(union(ens1,ens2))
#       }
#     }
#     
#     c1.list = apply(chevauche_1_2,1,which.max)
#     c2.list = apply(chevauche_1_2,2,which.max)
#     # Composantes du réplicat 1 et leur composante chevauchant le plus dans le réplicat 2
#     chev_crh_comp1 = rbind((1:structure_net_comp1[[b]]$no)[structure_net_comp1[[b]]$csize>1], unlist(c1.list[structure_net_comp1[[b]]$csize>1]))
#     # Composantes du réplicat 2 et leur composante chevauchant le plus dans le réplicat 1
#     chev_crh_comp2 = rbind((1:structure_net_comp2[[b]]$no)[structure_net_comp2[[b]]$csize>1], unlist(c2.list[structure_net_comp2[[b]]$csize>1]))
#     
#     
#     chev_crh_comp = matrix(NA, 2, ncol = min(ncol(chev_crh_comp1), ncol(chev_crh_comp2)))
#     
#     for (m1 in seq_len(ncol(chev_crh_comp1))) {
#       for (m2 in seq_len(ncol(chev_crh_comp2))) {
#         if(all(chev_crh_comp1[, m1]==rev(chev_crh_comp2[, m2]))){
#           chev_crh_comp[, sum(!is.na(chev_crh_comp[1, ])) + 1] <- chev_crh_comp1[, m1]
#         } 
#       }
#     }
#     
#     chev_crh_comp <- matrix(chev_crh_comp[!is.na(chev_crh_comp)], 2, ncol = length(chev_crh_comp[!is.na(chev_crh_comp)])/2)
#     
#     results[[length(results)+1]] <- list(
#       "chevauche_1_2_perc" = chevauche_1_2,
#       "chev_crh_comp1" = chev_crh_comp1,
#       "chev_crh_comp2" = chev_crh_comp2,
#       "chev_crh_comp" = chev_crh_comp,
#       "no" = c(length(unique(c1.list)), length(unique(c2.list)))
#     )
#     
#   }
#   
#   names(results) <- str_c("chr_overlap", seq_len(n_block))
#   results
# 
# }



