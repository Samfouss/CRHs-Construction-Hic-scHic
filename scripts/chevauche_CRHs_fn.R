
#' Cette fonction permet d'appréhender le taux de chevauchement entre les CRHs
#'
#' @description
#' #' `dist_func` retourne la distance entre deux points
#'
#' @param point_1 
#' @param point_2 

chr_overlap <- function(structure_net_comp1, structure_net_comp2){
  
  chevauche_1_2 = matrix(NA,structure_net_comp1$no,structure_net_comp2$no)
  for (i in (1:structure_net_comp1$no)[structure_net_comp1$csize>1])
  {
    for (j in (1:structure_net_comp2$no)[structure_net_comp2$csize>1])
    {
      ens1 = names(structure_net_comp1$membership)[structure_net_comp1$membership==i]
      ens2 = names(structure_net_comp2$membership)[structure_net_comp2$membership==j]
      chevauche_1_2[i,j] = length(intersect(ens1,ens2)) /length(union(ens1,ens2))
    }
  }
  
  c1.list = apply(chevauche_1_2,1,which.max)
  c2.list = apply(chevauche_1_2,2,which.max)
  # Composantes du réplicat 1 et leur composante chevauchant le plus dans le réplicat 2
  chev_chr_comp1 = rbind((1:structure_net_comp1$no)[structure_net_comp1$csize>1],
        unlist(c1.list[structure_net_comp1$csize>1]))
  # Composantes du réplicat 2 et leur composante chevauchant le plus dans le réplicat 1
  chev_chr_comp2 = rbind((1:structure_net_comp2$no)[structure_net_comp2$csize>1],
        unlist(c2.list[structure_net_comp2$csize>1]))
  
  results = list(
    "chevauche_1_2_perc" = chevauche_1_2,
    "chev_chr_comp1" = chev_chr_comp1,
    "chev_chr_comp2" = chev_chr_comp2
  )
}

