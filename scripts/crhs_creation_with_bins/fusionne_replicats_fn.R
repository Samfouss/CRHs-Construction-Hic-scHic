
###################### Function 1 ###################

#' @description Fusionne deux CRHs appariés dans des réplicats différents
#' @param structure_net_comp1 Composantes du réplicat 1 retournées par la fonction components
#' @param structure_net_comp2 Composantes du réplicat 2 retournées par la fonction components
#' @param dist_bin1 Matrice d'adjacence du réplicat 1
#' @param dist_bin2 Matrice d'adjacence du réplicat 2
#' @param edge_overlap_result On donne à la fonction fusion_comp les résulats retournés par la fonction edgeoverlap


bip_graphs_fusion <- function(structure_net_comp1, structure_net_comp2, edge_overlap_result){
  
  # On recupère les blocs à parcourir
  # n_block = max(
  #   length(structure_net_comp1), 
  #   length(structure_net_comp2)
  # )
  n_block = 16
  # Initialisation de la liste à recevoir les résultats
  results = list()
  
  for (b in seq_len(n_block)){
    
    if(b == 1){
      
      results[[length(results)+1]] <- list(
        "more_info" = list(),
        "crhs" = list()
      )
      
      next
    }
    
    crhs = list()
    nb_crh1 = length(structure_net_comp1[[b]]$crhs)
    nb_crh2 = length(structure_net_comp2[[b]]$crhs)
    
    paires.crh <- edge_overlap_result[[b]]$chev_edge_comp
    # print(paires.crh)
    # On boucle sur les CRHs à apparier
    
    rest_add1 = seq_len(nb_crh1)
    rest_add2 = seq_len(nb_crh2)
    
    # print(paires.crh)
    
    if(sum(paires.crh) != 0){
      for (j in seq_len(ncol(paires.crh))) {
        
        # On récupère ici les lignes des éléments qui appartiennent au CRHs à fusionner dans le premier replica
        ens1 = paires.crh[1, j]
        # On récupère ici les lignes des éléments qui appartiennent au CRHs à fusionner dans le deuxième replica
        ens2 = paires.crh[2, j]
        
        # Dans ens1 et ens2, on les noms des différents éléments des différents CRHs.
        mat1 = structure_net_comp1[[b]]$crhs[[ens1]]$mat_incidence
        mat2 = structure_net_comp2[[b]]$crhs[[ens2]]$mat_incidence
        
        rowmat = union(rownames(mat1), rownames(mat2))
        colmat = union(colnames(mat1), colnames(mat2))
        
        mat_merge_1 <- matrix(
          0,
          ncol=length(colmat), 
          nrow=length(rowmat), 
          dimnames=list(rowmat, colmat)
        )
        mat_merge_2 <- mat_merge_1
        
        indxA <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat1), colnames(mat1), FUN=paste)
        indxB <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat2), colnames(mat2), FUN=paste)
        mat_merge_1[indxA] <- mat1
        mat_merge_2[indxB] <- mat2
        
        # Ici les matrices mat_merge_1 et mat_merge_2 on les meme dimenssions avec les memes noms. Elles correspondent respectivement et exactement aux matrices adj_fusion et cc.list[[i]], juste qu'elles ont des dimenssions égales
        # On peut donc à ce stade faire la somme des deux matrices
        adj_fusion = mat_merge_1 + mat_merge_2
        # Une fois la somme faite, il faut ramener les valeurs qui sont supérieures ou égales à 2 (seulement dans notre cas d'étude la valeur maximale est 2) à 1.
        adj_fusion[adj_fusion>1] <- 1
        
        name1 = structure_net_comp1[[b]]$crhs[[ens1]]$name
        name2 = structure_net_comp2[[b]]$crhs[[ens2]]$name
        
        crhs[[length(crhs)+1]] = list(
          "name" = str_c(name1, " - ", name2),
          "mat_incidence" = adj_fusion
        )
        
      }
      # Une fois les CRHs fusionnés dans la boucle précédente, on doit juste à ce niveau ajouter ceux qui ne sont pas fusionnés avec leurs noms par la meme occasion
      rest_add1 = rest_add1[!rest_add1 %in% paires.crh[1, ]]
      rest_add2 = rest_add2[!rest_add2 %in% paires.crh[2, ]]
    }
    
    if(length(rest_add1) != 0){
      for (i in rest_add1) {
        crhs[[length(crhs)+1]] = list(
          "name" = structure_net_comp1[[b]]$crhs[[i]]$name,
          "mat_incidence" = structure_net_comp1[[b]]$crhs[[i]]$mat_incidence
        )
      }
    }
    
    
    if(length(rest_add2) != 0){
      for (j in rest_add2) {
        crhs[[length(crhs)+1]] = list(
          "name" = structure_net_comp2[[b]]$crhs[[j]]$name,
          "mat_incidence" = structure_net_comp2[[b]]$crhs[[j]]$mat_incidence
        )
      }
    }
    
    
    #################################################################"
    
    if(length(crhs)>0){
      names(crhs) <- str_c("crh", 1:length(crhs))
      
      # On a fini avec les CRHs individuellement, maintenant ce qui serait interessant de reconstituer toute la matrice d'incidence à partir des différents CRHs
      
      # - On commence par prendre la toute première matrice d'incidence
      # all_mat_incidence = crhs[[1]]$mat_incidence
      resume = matrix(NA, nrow = length(crhs), ncol = 1)
      resume[1, 1] <- crhs[[1]]$name
      
      if(length(crhs)>1){
        for (k in 2:length(crhs)) {
          
          resume[k, 1] <- crhs[[k]]$name
          
        } 
      }
      
      results[[length(results)+1]] <- list(
        "more_info" = resume,
        "crhs" = crhs
      )
    }
    
    # if(sum(structure_net_comp1[[b]]$crhs[[1]]$mat_incidence)==-1){
    #   crhs = structure_net_comp2[[b]]$crhs
    # }
    # 
    # if(sum(structure_net_comp2[[b]]$crhs[[1]]$mat_incidence)==-1){
    #   crhs = structure_net_comp1[[b]]$crhs
    # }
    
  }
  
  names(results) <- str_c("block", seq_len(length(results)))
  results
  
}


# fusion_comp2(net1, net2, res)
