#' Cette fonction permet d'appréhender le taux de chevauchement entre les CRHs
#'
#' @description
#' #' `dist_func` retourne la distance entre deux points
#'
#' @param point_1 
#' @param point_2 

chr_overlap <- function(structure_net_comp1, structure_net_comp2){
  
  n_block = length(structure_net_comp1) - 1
  results = list()
  
  for (b in seq_len(n_block)) {
    
    chevauche_1_2 = matrix(NA,structure_net_comp1[[b]]$no,structure_net_comp2[[b]]$no)
    for (i in (1:structure_net_comp1[[b]]$no)[structure_net_comp1[[b]]$csize>1]){
      for (j in (1:structure_net_comp2[[b]]$no)[structure_net_comp2[[b]]$csize>1]){
        ens1 = names(structure_net_comp1[[b]]$membership)[structure_net_comp1[[b]]$membership==i]
        ens2 = names(structure_net_comp2[[b]]$membership)[structure_net_comp2[[b]]$membership==j]
        chevauche_1_2[i,j] = length(intersect(ens1,ens2)) /length(union(ens1,ens2))
      }
    }
    
    c1.list = apply(chevauche_1_2,1,which.max)
    c2.list = apply(chevauche_1_2,2,which.max)
    # Composantes du réplicat 1 et leur composante chevauchant le plus dans le réplicat 2
    chev_crh_comp1 = rbind((1:structure_net_comp1[[b]]$no)[structure_net_comp1[[b]]$csize>1], unlist(c1.list[structure_net_comp1[[b]]$csize>1]))
    # Composantes du réplicat 2 et leur composante chevauchant le plus dans le réplicat 1
    chev_crh_comp2 = rbind((1:structure_net_comp2[[b]]$no)[structure_net_comp2[[b]]$csize>1], unlist(c2.list[structure_net_comp2[[b]]$csize>1]))
    
    
    chev_crh_comp = matrix(NA, 2, ncol = min(ncol(chev_crh_comp1), ncol(chev_crh_comp2)))
    
    for (m1 in seq_len(ncol(chev_crh_comp1))) {
      for (m2 in seq_len(ncol(chev_crh_comp2))) {
        if(all(chev_crh_comp1[, m1]==rev(chev_crh_comp2[, m2]))){
          chev_crh_comp[, sum(!is.na(chev_crh_comp[1, ])) + 1] <- chev_crh_comp1[, m1]
        } 
      }
    }
    
    chev_crh_comp <- matrix(chev_crh_comp[!is.na(chev_crh_comp)], 2, ncol = length(chev_crh_comp[!is.na(chev_crh_comp)])/2)
    
    results[[length(results)+1]] <- list(
      "chevauche_1_2_perc" = chevauche_1_2,
      "chev_crh_comp1" = chev_crh_comp1,
      "chev_crh_comp2" = chev_crh_comp2,
      "chev_crh_comp" = chev_crh_comp,
      "no" = c(length(unique(c1.list)), length(unique(c2.list)))
    )
    
  }
  
  names(results) <- str_c("chr_overlap", seq_len(n_block))
  results

}


#' Cette fonction permet d'appréhender le taux de chevauchement entre les CRHs
#'
#' @description
#' #' `edge_identity_overlap` retourne la similarité des CRHs en terme de nombre d'arrêts. L'idée est de partir de la matrice binaire des distance entre les différentes billes et de compter exactement le nombre de billes qui concident en terme de nombre d'arrêtes (correspond a 1 dans la matrice binaire) sur le nombre de billes total
#'
#' @param point_1 
#' @param point_2 
edge_identity_overlap <- function(structure_net_comp1, structure_net_comp2){
  
  n_block = length(structure_net_comp1) - 1
  results = list()
  
  for (b in seq_len(n_block)){
    
    dist_bin1 = structure_net_comp1[[b]]$dist_bin_bip
    dist_bin2 = structure_net_comp2[[b]]$dist_bin_bip
    chevauche_1_2 = matrix(
      NA,
      structure_net_comp1[[b]]$no_bip,
      structure_net_comp2[[b]]$no_bip
    )
    
    for (i in (1:structure_net_comp1[[b]]$no_bip)[structure_net_comp1[[b]]$csize_bip>1]){
      
      mat1 <- dist_bin1
      if(is.matrix(mat1)){
        mat1[mat1 != 0]<- 0
        
        membership1 <- structure_net_comp1[[b]]$membership_bip
        mat1[
          row.names(dist_bin1) %in% names(membership1[membership1==i]),
          names(membership1[membership1==i])
        ] = dist_bin1[
          row.names(dist_bin1) %in% names(membership1[membership1==i]),
          names(membership1[membership1==i])
        ]
      }
      
      for (j in (1:structure_net_comp2[[b]]$no_bip)[structure_net_comp2[[b]]$csize_bip>1]){
        
        mat2 <- dist_bin2
        if(is.matrix(mat2)){
          mat2[mat2 != 0]<- 0
          
          membership2 <- structure_net_comp2[[b]]$membership_bip
          mat2[
            row.names(dist_bin2) %in% names(membership2[membership2==i]),
            names(membership2[membership2==i])
          ] = dist_bin2[
            row.names(dist_bin2) %in% names(membership2[membership2==i]),
            names(membership2[membership2==i])
          ]
        }
        chevauche_1_2[i,j] = sum(mat1&mat2)/sum(mat1|mat2)
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


