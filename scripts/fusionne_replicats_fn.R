
fusion_comp <- function(structure_net_comp1, structure_net_comp2, edge_overlap_result){
  
  #' @description Fusionne deux CRHs appariés dans des réplicats différents
  #' @param structure_net_comp1 Composantes du réplicat 1 retournées par la fonction components
  #' @param structure_net_comp2 Composantes du réplicat 2 retournées par la fonction components
  #' @param dist_bin1 Matrice d'adjacence du réplicat 1
  #' @param dist_bin2 Matrice d'adjacence du réplicat 2
  #' @param edge_overlap_result On donne à la fonction fusion_comp les résulats retournés par la fonction edgeoverlap
  
  results = list()
  n_block <- length(structure_net_comp1)

  for (b in seq_len(n_block)) {
    
    cc.list <- list()
    dist_bin1 <- structure_net_comp1[[b]]$dist_bin
    dist_bin2 <- structure_net_comp2[[b]]$dist_bin
    paires.crh <- edge_overlap_result[[b]]$chev_edge_comp
    # Dans ces variables permettrons de voir si pour un block donnée dans un réplicat, il existe des promoters connectés à un bin au moins
    promoters_exist1 = sum(structure_net_comp1[[b]]$dist_bin)
    promoters_exist2 = sum(structure_net_comp2[[b]]$dist_bin)
    # On boucle sur les CRHs à apparier
    
    # Si on rencontre un block sans connection entre promoter et bins dans un replicat, on continue avec le second réplicat
    if (promoters_exist1==0){
      
      cc.list[[length(cc.list)+1]] <- dist_bin2
      #next
    }else if(promoters_exist2 == 0){
      
      cc.list[[length(cc.list)+1]] <- dist_bin1
      #next
    
    # Si les deux blocks des deux replicat n'ont pas de promoters, rien a faire. On initialise les valeur de cc.list à zéro afin de les utiliser après
    }else if(promoters_exist1 == 0 & promoters_exist2 == 0){
      
      cc.list[[length(cc.list)+1]] <- 0
      
    }else{
      for (j in seq_len(ncol(paires.crh))) {
        
        # On récupère ici les lignes des éléments qui appartiennent au CRHs à fusionner dans le premier replica
        ens1 = names(structure_net_comp1[[b]]$membership_bip)[structure_net_comp1[[b]]$membership_bip==paires.crh[1, j]]
        # On récupère ici les lignes des éléments qui appartiennent au CRHs à fusionner dans le deuxième replica
        ens2 = names(structure_net_comp2[[b]]$membership_bip)[structure_net_comp2[[b]]$membership_bip==paires.crh[2, j]]
        noms = union(ens1,ens2)
        # On crée une matrice où on doit mettre les valeurs de matrice bipartite des éléments ens1 et ens2. Cette matrice n'est rien d'autre que l'union des deux ensembles
        adjcomp = matrix(0,length(noms),length(noms))
        dimnames(adjcomp) = list(noms,noms)
        
        adjcomp[ens1,ens1] = dist_bin1[ens1,ens1]
        adjcomp[ens2,ens2] = adjcomp[ens2,ens2]|dist_bin2[ens2,ens2]
        
        cc.list[[length(cc.list)+1]] <- adjcomp
        
      }
    }

    # names(cc.list) <- str_c("matched", seq_len(ncol(paires.crh)))
    # 
    # results[[length(results)+1]] <- cc.list
    
    # On ajoute les composantes non appariées
    j = length(cc.list)
    
    if (promoters_exist1 != 0 & promoters_exist2 != 0){
      # du réplicat 1
      for(i in (1:structure_net_comp1[[b]]$no_bip)[!(1:structure_net_comp1[[b]]$no_bip)%in%paires.crh[1, ] & structure_net_comp1[[b]]$csize_bip>1]){
        j = j+1
        cc.list[[j]] = structure_net_comp1[[b]]$dist_bin[
          names(structure_net_comp1[[b]]$membership_bip==i),
          names(structure_net_comp1[[b]]$membership_bip==i)
        ]
      }
      # du réplicat 2
      for(i in (1:structure_net_comp2[[b]]$no_bip)[!(1:structure_net_comp2[[b]]$no_bip)%in%paires.crh[2, ] & structure_net_comp2[[b]]$csize_bip>1]){
        j = j+1
        cc.list[[j]] = structure_net_comp2[[b]]$dist_bin[
          names(structure_net_comp2[[b]]$membership_bip==i),
          names(structure_net_comp2[[b]]$membership_bip==i)
        ]
      }
      
      ## Création d'une nouvelle matrice d'adjacence
      # nr.vec = c(0,cumsum(sapply(cc.list,nrow)))
      # nr = nr.vec[length(nr.vec)]
      # adj_fusion = matrix(0,nr,nr)
      # row_name = c()
      # col_name = c()
      # for (i in 1:length(cc.list)){
      #   adj_fusion[
      #     (nr.vec[i]+1):nr.vec[i+1],
      #     (nr.vec[i]+1):nr.vec[i+1]
      #   ] = cc.list[[i]]
      #   row_name <- c(row_name, rownames(cc.list[[i]]))
      #   col_name <- c(col_name, colnames(cc.list[[i]]))
      # }
      # rownames(adj_fusion) <- row_name
      # colnames(adj_fusion) <- col_name


      # Création d'une nouvelle matrice d'adjacence à partir des noms des matrices, histoire qu'il n'y ait pas de répétitions
      # 1 - On initialise la matrice adj_fusion avec la première matrice de cc.list
      adj_fusion = cc.list[[1]]
      
      if (length(cc.list)>1){
        for (i in 2:length(cc.list)){
          # On initialise les matrice à fusionner dans mat1 et mat2
          mat1 = adj_fusion
          mat2 = cc.list[[i]]
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
          
        }
      }
    }
    
    net_bip <- graph_from_incidence_matrix(
      adj_fusion
    )
    net_components_bip <- components(net_bip, mode = c("weak", "strong"))
    
    plot_main_cluster <- which(net_components_bip$csize>1)
    vert_ids <- V(net_bip)[net_components_bip$membership %in% plot_main_cluster]
    net_to_plot <- induced_subgraph(net_bip, vert_ids)
    V(net_to_plot)$color <- V(net_to_plot)$type
    V(net_to_plot)$color=gsub("FALSE","red",V(net_to_plot)$color)
    V(net_to_plot)$color=gsub("TRUE","lightblue",V(net_to_plot)$color)
    
    file_name = paste0("crhs_merge_block_", b, ".jpeg")
    jpeg(file_name, width = 980, height = 680)
    par(mar = c(1, 1, 1, 1))
    plot(net_to_plot, edge.arrow.size=.2,vertex.label=NA, main = paste0("Représenation graphique du block ", b, " des CRHs fusionnés "))
    dev.off()
    
    results[[length(results)+1]] <- list(
      "dist_bin" = adj_fusion,
      "net_to_plot" = net_to_plot,
      "membership_bip" = net_components_bip$membership,
      "csize_bip" = net_components_bip$csize,
      "no_bip" = net_components_bip$no
    ) 
  }
  
  names(results) <- str_c("block", seq_len(n_block))
  results
  
}




#' fusion_comp <- function(structure_net_comp1, structure_net_comp2, edge_overlap_result){
#'   
#'   #' @description Fusionne deux CRHs appariés dans des réplicats différents
#'   #' @param structure_net_comp1 Composantes du réplicat 1 retournées par la fonction components
#'   #' @param structure_net_comp2 Composantes du réplicat 2 retournées par la fonction components
#'   #' @param edge_overlap_result Vecteur de longueur 2 avec l'indice de la composante du réplicat 1 et de celle du réplicat 2 à fusionner
#'   
#'   results = list()
#'   n_block <- length(structure_net_comp1)
#' 
#'   for (b in seq_len(n_block)) {
#'     
#'     paires.crh <- edge_overlap_result[[b]]$chev_edge_comp
#'     
#'     for (p in seq_len(ncol(paires.crh))) {
#'       
#'       # Matrice d'adjacence du réplicat 1
#'       dist_bin1 <- structure_net_comp1[[b]]$dist_bin
#'       # Matrice d'adjacence du réplicat 2
#'       dist_bin2 <- structure_net_comp2[[b]]$dist_bin
#'       
#'       ens1 = names(structure_net_comp1[[b]]$membership_bip)[structure_net_comp1[[b]]$membership_bip==paires.crh[1, p]]
#'       ens2 = names(structure_net_comp2[[b]]$membership_bip)[structure_net_comp2[[b]]$membership_bip==paires.crh[2, p]]
#'       noms = union(ens1,ens2)
#'       adjcomp = matrix(0,length(noms),length(noms))
#'       dimnames(adjcomp) = list(noms,noms)
#'       
#'       adjcomp[ens1,ens1] = dist_bin1[ens1,ens1]
#'       adjcomp[ens2,ens2] = adjcomp[ens2,ens2]|dist_bin2[ens2,ens2]
#'       
#'       results[[length(results)+1]] <- list(
#'         "adjcomp" = adjcomp,
#'         "ens1" = ens1,
#'         "ens2" = ens2
#'       )
#'       
#'     }
#'   }
#'   
#'   names(results) <- str_c("adjcomp", seq_len(n_block))
#'   results
#'   
#' }


