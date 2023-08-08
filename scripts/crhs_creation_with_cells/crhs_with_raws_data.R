
library(igraph)

# Chargement des promoters
load("rdata/all_rda_data/scHic_promoters_ids.rda")

ncells = 250
mat_dim = 562

raw_data_crhs = list()

for (cell in seq_len(ncells)) {
  
  cell_matrix = read.table(paste0("rdata/single_cell_hic_data/hic_ideal_data_with_rep/hic_mat_cell_", sprintf("%03d", cell), ".txt"))
  
  cell_matrix <- sapply(
    cell_matrix, 
    FUN = function(item) ifelse(item>=1, 1, 0)
  )
  name = paste0("BAC", sprintf("%03d", seq_len(mat_dim)))
  row.names(cell_matrix) = name
  colnames(cell_matrix) = name
  
  cell_matrix <- cell_matrix[
    scHic_promoters_ids, 
    -scHic_promoters_ids, 
    drop=FALSE
  ]
  
  
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

  crhs = list()
  for (i in (1:net_components_bip$no)[net_components_bip$csize>1]){
    
    members <- net_components_bip$membership
    mat_bin = cell_matrix[
      rownames(cell_matrix) %in% names(members[members==i]),
      colnames(cell_matrix) %in% names(members[members==i]),
      drop=FALSE
    ]
    
    crhs[[length(crhs)+1]] = list(
      "name" = paste0("( ", cell, " )"),
      "mat_incidence" = mat_bin
    )
    
  }
  
  names(crhs) <- paste0("crh", 1:length(crhs))
  
  res <- list(
    "dist_bin" = cell_matrix,
    "net_to_plot" = net_to_plot,
    "membership_bip" = net_components_bip$membership,
    "csize_bip" = net_components_bip$csize,
    "no_bip" = net_components_bip$no,
    "crhs" = crhs
  )
  
  raw_data_crhs[[length(raw_data_crhs)+1]] <- res
  
  names(raw_data_crhs) <- paste0("raw_data", seq_len(length(raw_data_crhs)))
}



nb_crhs = 0
for (cell in seq_len(ncells)) {
  for (i in seq_len(length(raw_data_crhs[[cell]]$crhs))) {
    nb_crhs = nb_crhs + length(raw_data_crhs[[cell]]$crhs[[i]])
  }
}
nb_crhs

crhs_inspection = matrix(
  0,
  nrow = nb_crhs,
  ncol = 3
)

l = 1

for (cell in seq_len(ncells)) {
  for (i in seq_len(length(raw_data_crhs[[cell]]$crhs))) {
    crhs_inspection[l, 1] =  dim(raw_data_crhs[[cell]]$crhs[[i]]$mat_incidence)[1]
    crhs_inspection[l, 2] =  dim(raw_data_crhs[[cell]]$crhs[[i]]$mat_incidence)[2]
    crhs_inspection[l, 3] =  paste0("cell ", sprintf("%03d", cell), "- CRH ", i)
    l = l +  1
  }
}

table(crhs_inspection[])

