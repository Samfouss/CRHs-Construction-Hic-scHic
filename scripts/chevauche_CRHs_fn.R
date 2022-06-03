
#' Cette fonction permet d'appréhender le taux de chevauchement entre les CRHs
#'
#' @description
#' #' `dist_func` retourne la distance entre deux points
#'
#' @param point_1 
#' @param point_2 

chr_overlap <- function(structure_net_comp1, structure_net_comp2){
  
  chevauche_1_2 = matrix(NA,structure_net_comp1$no,structure_net_comp2$no)
  for (i in (1:structure_net_comp1$no)[structure_net_comp1$csize>1]){
    for (j in (1:structure_net_comp2$no)[structure_net_comp2$csize>1]){
      ens1 = names(structure_net_comp1$membership)[structure_net_comp1$membership==i]
      ens2 = names(structure_net_comp2$membership)[structure_net_comp2$membership==j]
      chevauche_1_2[i,j] = length(intersect(ens1,ens2)) /length(union(ens1,ens2))
    }
  }
  
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
}

edge_identity_overlap <- function(structure_net_comp1, structure_net_comp2, dist_bin1, dist_bin2){
  
  chevauche_1_2 = matrix(NA,structure_net_comp1$no,structure_net_comp2$no)
  for (i in (1:structure_net_comp1$no)[structure_net_comp1$csize>1]){
    mat1 = matrix(0,nrow(dist_bin1),ncol(dist_bin1))
    mat1[structure_net_comp1$membership==i,structure_net_comp1$membership==i] = dist_bin1[structure_net_comp1$membership==i,structure_net_comp1$membership==i]
    for (j in (1:structure_net_comp2$no)[structure_net_comp2$csize>1]){
      mat2 = matrix(0,nrow(dist_bin2),ncol(dist_bin2))
      mat2[structure_net_comp1$membership==j,structure_net_comp1$membership==j] = dist_bin2[structure_net_comp1$membership==j,structure_net_comp1$membership==j]
      chevauche_1_2[i,j] = sum(mat1&mat2) /sum(mat1|mat2)
    }
  }
  
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
}


edges_overlap <- function(data_frame1, data_frame2, rep1, rep2, edge_var, crh_var){
  
  crh_list1 <- unique(data_frame1[, crh_var])
  crh_list2 <- unique(data_frame2[, crh_var])
  
  chevauche = matrix(NA,length(unique(data_frame1[, crh_var])),length(unique(data_frame2[, crh_var])))
  
  rownames(chevauche) <- paste0("R", rep1, "B", unique(data_frame1[, "X4"]), "CRH", unique(data_frame1[, crh_var]))
  colnames(chevauche) <- paste0("R", rep2, "B", unique(data_frame2[, "X4"]), "CRH", unique(data_frame2[, crh_var]))
  
  
  for (i in seq_len(length(crh_list1))) {
    ens1 <- filter(data_frame1, eval(parse(text=crh_var))==crh_list1[i])[, edge_var]
    for (j in seq_len(length(crh_list2))) {
      ens2 <- filter(data_frame2, eval(parse(text=crh_var))==crh_list2[j])[, edge_var]
      chevauche[i,j] <- length(intersect(ens1,ens2)) /length(union(ens1,ens2))
    }
  }
  
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
}








