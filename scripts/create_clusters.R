
#' Cette fonction permet d'appréhender le taux de chevauchement entre les CRHs
#'
#' @description
#' #' `create_clust_graph` 
#'
#' @param all_data 
#' @param rep_num 


create_clust_graph <- function(all_data, rep_num, block, min_dist = 3){
  
  par_default <- par(bty = 'n')
  
  three_dim_data <- all_data[all_data$X4==block, 1:3]
  row.names(three_dim_data) <- pull(all_data[all_data$X4==block, 5])
  
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
      pull(all_data[all_data$X4==block, 5]),
      pull(all_data[all_data$X4==block, 5])
    )
  )
  
  net <- graph_from_adjacency_matrix(
    compute_dist_bin, 
    mode='undirected'
  )
  net <- simplify(net, remove.multiple = FALSE, remove.loops = TRUE)
  # dev.off()
  par(mar = c(1, 1, 1, 1))
  g <- plot(net, edge.arrow.size=.4,vertex.label=NA, main = paste0("Représenation graphique du block ", block, " pour le replicat ", rep_num))
  g
  par(par_default)
  
  edges <- data.frame(table(unlist(strsplit(attr(E(net)[which_mutual(net)],"vnames"),"\\|"))))%>%
    rename(
      ID = Var1,
      nb_edges = Freq
    )
  
  # , vertex.color="yellow",vertex.size=15,vertex.label.color="blue",edge.color="red",edge.width=5
  
  res = list(
    #"dist" = compute_dist,
    #"dist_bin" = compute_dist_bin,
    "edges" = edges,
    "plot" = g,
    "net" = net
  )
}
