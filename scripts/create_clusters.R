


################# function 2 ##################
create_bip_clust_graph2 <- function(all_data, promoters_vec, rep_num, block_vec, min_dist = 3){
  
  par_default <- par(bty = 'n')
  res = list()
  block <- unique(block_vec)
  
  # structure = data.frame()
  
  for(i in 1:length(block)) {

    b = block[i]
    all_data <- all_data%>%
      mutate(
        is_promoter = if_else(row_number() %in% promoters_vec , 1, 0)
      )
    
    three_dim_data <- all_data[all_data$X4==b, 1:3]
    block_promoters <- which(
      all_data[all_data$X4==b, "is_promoter"]$is_promoter==1
    )
    
    # S'il n y a pas de promoteurs dans le block
    if(length(block_promoters)!=0){
      
      mat_row_name = pull(all_data[all_data$X4==b, 5])
      #three_dim_data = column_to_rownames(three_dim_data, var = "ID")
      
      compute_dist <- as.matrix(
        dist(
          x = three_dim_data, 
          method = "euclidean", 
          diag = TRUE
        )
      )
      
      # 342 billes
      # 13 lines  promoters (342 lignes : 342 - 13 lignes à 0)
      # 342 colonnes
      compute_dist_bin <- sapply(
        compute_dist, 
        FUN = function(item) ifelse(item<=min_dist, 1, 0)
      )
      
      compute_dist_bin <- matrix(
        compute_dist_bin, 
        nrow = nrow(compute_dist), 
        ncol = ncol(compute_dist), 
        byrow = TRUE,
        dimnames = list(mat_row_name, mat_row_name)
      )
      
      diag(compute_dist_bin) <- 0
      
      compute_dist_bin <- compute_dist_bin[block_promoters, -block_promoters, drop=FALSE]
      
      #compute_dist_bin[!(1:nrow(compute_dist_bin) %in% block_promoters), ] <- 0
      # net <- graph_from_incidence_matrix(
      #   compute_dist_bin
      # )
      
      net_bip <- graph_from_incidence_matrix(
        compute_dist_bin
      )
      
      net_bip <- simplify(net_bip, remove.multiple = FALSE, remove.loops = TRUE)
      # net <- simplify(net, remove.multiple = FALSE, remove.loops = TRUE)
      
      net_components_bip <- components(net_bip, mode = c("weak", "strong"))
      # net_components <- components(net, mode = c("weak", "strong"))
      #biggest_cluster <- which.max(net_components$csize)
      #net_comp2 <- clusters(net, mode="weak")
      plot_main_cluster <- which(net_components_bip$csize>1)
      vert_ids <- V(net_bip)[net_components_bip$membership %in% plot_main_cluster]
      net_to_plot <- induced_subgraph(net_bip, vert_ids)
      V(net_to_plot)$color <- V(net_to_plot)$type
      V(net_to_plot)$color=gsub("FALSE","red",V(net_to_plot)$color)
      V(net_to_plot)$color=gsub("TRUE","lightblue",V(net_to_plot)$color)
      
      # file_name = paste0("graphique_rep", rep_num, "_block", b, ".jpeg")
      # jpeg(file_name, width = 980, height = 680)
      # par(mar = c(1, 1, 1, 1))
      # plot(net_to_plot, edge.arrow.size=.2,vertex.label=NA, main = paste0("Représenation graphique du block ", b, " pour le replicat ", rep_num))
      # dev.off()
      
      
      
      # A ce niveau, on prend chaque matrice d'incidence des différents CRHs qu'on ajoute à une liste 
      crhs = list()
      for (i in (1:net_components_bip$no)[net_components_bip$csize>1]){
        
        members <- net_components_bip$membership
        mat_bin = compute_dist_bin[
          rownames(compute_dist_bin) %in% names(members[members==i]),
          colnames(compute_dist_bin) %in% names(members[members==i]),
          drop=FALSE
        ]
        
        crhs[[length(crhs)+1]] = list(
          "name" = str_c("( ", rep_num, ", ", i, " )"),
          "mat_incidence" = mat_bin
        )
        
      }
      
      # Il y a des blocks qui ont des promoteurs par contre aucune bille n'est connetées à ce promoteur
      
      if(length(crhs)>0){
        
        names(crhs) <- str_c("crh", 1:length(crhs))
        
        res[[length(res)+1]] <- list(
          "dist_bin" = compute_dist_bin,
          "net_to_plot" = net_to_plot,
          "membership_bip" = net_components_bip$membership,
          "csize_bip" = net_components_bip$csize,
          "no_bip" = net_components_bip$no,
          "crhs" = crhs
        )
        
      }else{
        
        crhs[[length(crhs)+1]] = list(
          "name" = str_c("( ", rep_num, ", ", i, " )"),
          "mat_incidence" = -1
        )
        
        names(crhs) <- str_c("crh", 1:length(crhs))
        
        res[[length(res)+1]] <- list(
          "dist_bin" = -1,
          "net_to_plot" = -1,
          "membership_bip" = -1,
          "csize_bip" = -1,
          "no_bip" = -1,
          "crhs" = crhs
        )
        
      }
      
      
    }else{
      
      crhs = list()
      crhs[[length(crhs)+1]] = list(
        "name" = str_c("( ", rep_num, ", ", i, " )"),
        "mat_incidence" = -1
      )
      
      names(crhs) <- str_c("crh", 1:length(crhs))
      
      res[[length(res)+1]] <- list(
        "dist_bin" = -1,
        "net_to_plot" = -1,
        "membership_bip" = -1,
        "csize_bip" = -1,
        "no_bip" = -1,
        "crhs" = crhs
      )
    }
    
  }
  
  # res[[length(res)+1]] <- list(
  #   "structure" = structure
  # )
  
  names(res) <- str_c("clust_block", seq_len(length(res)))
  res
  
}



