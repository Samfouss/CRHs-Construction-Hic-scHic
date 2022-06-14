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
  
<<<<<<< HEAD
  names(results) <- str_c("chr_overlap", seq_len(n_block))
  results
=======
  c1.list = apply(chevauche_1_2,1,which.max)
  c2.list = apply(chevauche_1_2,2,which.max)
  # Composantes du réplicat 1 et leur composante chevauchant le plus dans le réplicat 2
  chev_crh_comp1 = rbind((1:structure_net_comp1$no)[structure_net_comp1$csize>1], unlist(c1.list[structure_net_comp1$csize>1]))
  # Composantes du réplicat 2 et leur composante chevauchant le plus dans le réplicat 1
  chev_crh_comp2 = rbind((1:structure_net_comp2$no)[structure_net_comp2$csize>1], unlist(c2.list[structure_net_comp2$csize>1]))
  
  results = list(
    "chevauche_1_2_perc" = chevauche_1_2,
    "chev_crh_comp1" = chev_crh_comp1,
    "chev_crh_comp2" = chev_crh_comp2,
    "no" = c(length(unique(c1.list)), length(unique(c2.list)))
  )
>>>>>>> f17f9cda862cba86b7b0c7ba5053125c8850617a
}

edge_identity_overlap <- function(structure_net_comp1, structure_net_comp2, dist_bin1, dist_bin2){
  
  chevauche_1_2 = matrix(NA,structure_net_comp1$no,structure_net_comp2$no)
  for (i in (1:structure_net_comp1$no)[structure_net_comp1$csize>1]){
    mat1 = matrix(0,nrow(dist_bin1),ncol(dist_bin1))
    mat1[structure_net_comp1$membership==i,structure_net_comp1$membership==i] = dist_bin1[structure_net_comp1$membership==i,structure_net_comp1$membership==i]
    for (j in (1:structure_net_comp2$no)[structure_net_comp2$csize>1]){
      mat2 = matrix(0,nrow(dist_bin2),ncol(dist_bin2))
      mat2[structure_net_comp2$membership==j,structure_net_comp2$membership==j] = dist_bin2[structure_net_comp2$membership==j,structure_net_comp2$membership==j]
      chevauche_1_2[i,j] = sum(mat1&mat2) /sum(mat1|mat2)
    }
  }
  
  c1.list = apply(chevauche_1_2,1,which.max)
  c2.list = apply(chevauche_1_2,2,which.max)
  # Composantes du réplicat 1 et leur composante chevauchant le plus dans le réplicat 2
  chev_crh_comp1 = cbind((1:structure_net_comp1$no)[structure_net_comp1$csize>1], unlist(c1.list[structure_net_comp1$csize>1]))
  # Composantes du réplicat 2 et leur composante chevauchant le plus dans le réplicat 1
  chev_crh_comp2 = cbind((1:structure_net_comp2$no)[structure_net_comp2$csize>1], unlist(c2.list[structure_net_comp2$csize>1]))
  
  results = list(
    "chevauche_1_2_perc" = chevauche_1_2,
    "chev_crh_comp1" = chev_crh_comp1,
    "chev_crh_comp2" = chev_crh_comp2,
    "no" = c(length(unique(c1.list)), length(unique(c2.list)))
  )
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
    
    dist_bin1 = structure_net_comp1[[b]]$dist_bin
    dist_bin2 = structure_net_comp2[[b]]$dist_bin
    chevauche_1_2 = matrix(NA,structure_net_comp1[[b]]$no,structure_net_comp2[[b]]$no)
    
    for (i in (1:structure_net_comp1[[b]]$no)[structure_net_comp1[[b]]$csize>1]){
      mat1 = matrix(0,nrow(dist_bin1),ncol(dist_bin1))
      mat1[
        structure_net_comp1[[b]]$membership==i,
        structure_net_comp1[[b]]$membership==i] = dist_bin1[
          structure_net_comp1[[b]]$membership==i,
          structure_net_comp1[[b]]$membership==i
        ]
      
      for (j in (1:structure_net_comp2[[b]]$no)[structure_net_comp2[[b]]$csize>1]){
        mat2 = matrix(0,nrow(dist_bin2),ncol(dist_bin2))
        mat2[
          structure_net_comp2[[b]]$membership==j,
          structure_net_comp2[[b]]$membership==j] = dist_bin2[
            structure_net_comp2[[b]]$membership==j,
            structure_net_comp2[[b]]$membership==j
          ]
        chevauche_1_2[i,j] = sum(mat1&mat2)/sum(mat1|mat2)
      }
      
    }
    
    c1.list = apply(chevauche_1_2,1,which.max)
    c2.list = apply(chevauche_1_2,2,which.max)
    # Composantes du réplicat 1 et leur composante chevauchant le plus dans le réplicat 2
    chev_edge_comp1 = rbind(
      (1:structure_net_comp1[[b]]$no)[structure_net_comp1[[b]]$csize>1], 
      unlist(c1.list[structure_net_comp1[[b]]$csize>1])
    )
    # Composantes du réplicat 2 et leur composante chevauchant le plus dans le réplicat 1
    chev_edge_comp2 = rbind(
      (1:structure_net_comp2[[b]]$no)[structure_net_comp2[[b]]$csize>1], 
      unlist(c2.list[structure_net_comp2[[b]]$csize>1])
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
  
<<<<<<< HEAD
  names(results) <- str_c("edge_overlap", seq_len(n_block))
  results
  
=======
#  ens1_list = apply(chevauche,1,which.max)
#  ens2_list = apply(chevauche,2,which.max)
  ens1_list = dimnames(chevauche)[[2]][apply(chevauche, MARGIN = 1, FUN = which.max)]
  ens2_list = dimnames(chevauche)[[1]][apply(chevauche, MARGIN = 2, FUN = which.max)]
  
  # Composantes du réplicat 1 et leur composante chevauchant le plus dans le réplicat 2
  edges_chev_crh_comp1 = rbind(crh_list1, ens1_list)
  colnames(edges_chev_crh_comp1) <- rownames(chevauche)
  
  # Composantes du réplicat 2 et leur composante chevauchant le plus dans le réplicat 1
  edges_chev_crh_comp2 = rbind(crh_list2, ens2_list)
  colnames(edges_chev_crh_comp2) <- colnames(chevauche)
  
  results = list(
    "chevauche" = chevauche,
    "edges_chev_crh_comp1" = edges_chev_crh_comp1,
    "edges_chev_crh_comp2" = edges_chev_crh_comp2,
    "no" = c(length(unique(ens1_list)), length(unique(ens2_list)))
  )
>>>>>>> f17f9cda862cba86b7b0c7ba5053125c8850617a
}


