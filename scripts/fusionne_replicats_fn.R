<<<<<<< HEAD
fusion_comp <- function(structure_net_comp1, structure_net_comp2, dist_bin1, dist_bin2,paire){
=======
fusion_comp = function(structure_net_comp1, structure_net_comp2, dist_bin1, dist_bin2,paire){
>>>>>>> f17f9cda862cba86b7b0c7ba5053125c8850617a
  #' @description Fusionne deux CRHs appariés dans des réplicats différents
  #' @param structure_net_comp1 Composantes du réplicat 1 retournées par la fonction components
  #' @param structure_net_comp2 Composantes du réplicat 2 retournées par la fonction components
  #' @param dist_bin1 Matrice d'adjacence du réplicat 1
  #' @param dist_bin2 Matrice d'adjacence du réplicat 2
  #' @param paire Vecteur de longueur 2 avec l'indice de la composante du réplicat 1 et de celle du réplicat 2 à fusionner
<<<<<<< HEAD
  
=======
 
>>>>>>> f17f9cda862cba86b7b0c7ba5053125c8850617a
  ens1 = names(structure_net_comp1$membership)[structure_net_comp1$membership==paire[1]]
  ens2 = names(structure_net_comp2$membership)[structure_net_comp2$membership==paire[2]]
  noms = union(ens1,ens2)
  adjcomp = matrix(0,length(noms),length(noms))
  dimnames(adjcomp) = list(noms,noms)
<<<<<<< HEAD
  
=======

>>>>>>> f17f9cda862cba86b7b0c7ba5053125c8850617a
  adjcomp[ens1,ens1] = dist_bin1[ens1,ens1]
  adjcomp[ens2,ens2] = adjcomp[ens2,ens2]|dist_bin2[ens2,ens2]
  adjcomp
}

<<<<<<< HEAD

structure_1_net$clust_block1$membership[structure_1_net$clust_block1$membership==overlap$chr_overlap1$chev_crh_comp[1]]

=======
>>>>>>> f17f9cda862cba86b7b0c7ba5053125c8850617a
