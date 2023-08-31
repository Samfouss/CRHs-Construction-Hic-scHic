
call_igraph_fn <- function(cell_matrix, cell_num, resolution){
  
  net_bip <- graph_from_incidence_matrix(
    cell_matrix
  )
  
  # Effacer les connections multiples entre les noeuds
  if(!is_simple(net_bip)){
    net_bip <- simplify(net_bip, remove.multiple = FALSE, remove.loops = TRUE)
    # net <- simplify(net, remove.multiple = FALSE, remove.loops = TRUE) 
  }
  
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
    mat_bin = cell_matrix[
      rownames(cell_matrix) %in% names(members[members==i]),
      colnames(cell_matrix) %in% names(members[members==i]),
      drop=FALSE
    ]
    
    crhs[[length(crhs)+1]] = list(
      "name" = str_c("( ", cell_num, ".", resolution, " )"),
      "mat_incidence" = mat_bin
    )
    
  }
  
  names(crhs) <- str_c("crhs", ".", resolution, 1:length(crhs))

  res <- list(
    "dist_bin" = cell_matrix,
    "net_to_plot" = net_to_plot,
    "membership_bip" = net_components_bip$membership,
    "csize_bip" = net_components_bip$csize,
    "no_bip" = net_components_bip$no,
    "crhs" = crhs
  )
  
  res
  
}


