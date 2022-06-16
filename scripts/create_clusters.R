
#' Cette fonction permet de créer un igraph bipartite
#'
#' @description
#' #' `create_bip_clust_graph` 
#'
#' @param all_data 
#' @param rep_num 

create_bip_clust_graph <- function(all_data, promoters_vec, rep_num, block_vec, min_dist = 3){
  
  par_default <- par(bty = 'n')
  res = list()
  block <- unique(block_vec)
  
  structure = data.frame()
  
  for (i in 1:length(block)) {
    
    b = block[i]
    
    all_data <- all_data%>%
      mutate(
        is_promoter = if_else(row_number() %in% promoters_vec , 1, 0)
      )
    
    three_dim_data <- all_data[all_data$X4==b, 1:3]
    block_promoters <- which(all_data[all_data$X4==b, "is_promoter"]$is_promoter==1)
    row.names(three_dim_data) <- pull(all_data[all_data$X4==b, 5])
    
    compute_dist <- as.matrix(
      dist(
        x = three_dim_data, 
        method = "euclidean", 
        diag = TRUE
      )
    )
    
    compute_dist_bin <- sapply(
      compute_dist, 
      FUN = function(item) ifelse(item<=min_dist, 1, 0)
    )
    
    compute_dist_bin <- matrix(
      compute_dist_bin, 
      nrow = nrow(compute_dist), 
      ncol = ncol(compute_dist), 
      byrow = TRUE
    )
    
    compute_dist_bin_bip <- compute_dist_bin[block_promoters, ]
    
    net <- graph_from_incidence_matrix(
      compute_dist_bin_bip
    )
    net <- simplify(net, remove.multiple = FALSE, remove.loops = TRUE)
    # dev.off()
    par(mar = c(1, 1, 1, 1))
    g <- plot(net, edge.arrow.size=.4,vertex.label=NA, main = paste0("Représenation graphique du block ", b, " pour le replicat ", rep_num))
    g
    par(par_default)
    
    net_components <- components(net, mode = c("weak", "strong"))
    
    edges <- data.frame(table(unlist(strsplit(attr(E(net)[which_mutual(net)],"vnames"),"\\|"))))%>%
      rename(
        ID = Var1,
        nb_edges = Freq
      )
    
    structure <- rbind(
      structure, 
      edges%>%
        full_join(
          data.frame(
            ID = names(net_components$membership),
            crh_id = net_components$membership
          ),
          by = "ID"
        )
    )
    
    res[[length(res)+1]] <- list(
      "edges" = edges,
      "plot" = g,
      "dist_bin" = compute_dist_bin,
      "dist_bin_bip" = compute_dist_bin_bip,
      "membership" = net_components$membership,
      "csize" = net_components$csize,
      "no" = net_components$no
    )
    
  }
  
  res[[length(res)+1]] <- list(
    "structure" = structure
  )
  
  names(res) <- c(str_c("clust_block", block), "clust_block")
  res
  
}


#' Cette fonction permet d'appréhender le taux de chevauchement entre les CRHs
#'
#' @description
#' #' `create_clust_graph` 
#'
#' @param all_data 
#' @param rep_num 


create_clust_graph <- function(all_data, rep_num, block_vec, min_dist = 3){
  
  par_default <- par(bty = 'n')
  res = list()
  block <- unique(block_vec)
  
  structure = data.frame()
  
  for (i in 1:length(block)) {
    
    b = block[i]
    three_dim_data <- all_data[all_data$X4==b, 1:3]
    row.names(three_dim_data) <- pull(all_data[all_data$X4==b, 5])
    
    compute_dist <- as.matrix(
      dist(
        x = three_dim_data, 
        method = "euclidean", 
        diag = TRUE
      )
    )
    
    compute_dist_bin <- sapply(
      compute_dist, 
      FUN = function(item) ifelse(item<=min_dist, 1, 0)
    )

    
    compute_dist_bin <- matrix(
      compute_dist_bin, 
      nrow = nrow(compute_dist), 
      ncol = ncol(compute_dist), 
      byrow = TRUE,
      dimnames = list(
        pull(all_data[all_data$X4==b, 5]),
        pull(all_data[all_data$X4==b, 5])
      )
    )
    
    net <- graph_from_adjacency_matrix(
      compute_dist_bin, 
      mode='undirected'
    )
    net <- simplify(net, remove.multiple = FALSE, remove.loops = TRUE)
    # dev.off()
    par(mar = c(1, 1, 1, 1))
    g <- plot(net, edge.arrow.size=.4,vertex.label=NA, main = paste0("Représenation graphique du block ", b, " pour le replicat ", rep_num))
    g
    par(par_default)
    
    net_components <- components(net, mode = c("weak", "strong"))
    
    edges <- data.frame(table(unlist(strsplit(attr(E(net)[which_mutual(net)],"vnames"),"\\|"))))%>%
      rename(
        ID = Var1,
        nb_edges = Freq
      )
    
    structure <- rbind(
      structure, 
      edges%>%
        full_join(
          data.frame(
            ID = names(net_components$membership),
            crh_id = net_components$membership
          ),
          by = "ID"
        )
      )
    
    res[[length(res)+1]] <- list(
      "edges" = edges,
      "plot" = g,
      "dist_bin" = compute_dist_bin,
      "membership" = net_components$membership,
      "csize" = net_components$csize,
      "no" = net_components$no
    )
      
  }
  
  res[[length(res)+1]] <- list(
    "structure" = structure
  )
  
  names(res) <- c(str_c("clust_block", block), "clust_block")
  res
  
}

