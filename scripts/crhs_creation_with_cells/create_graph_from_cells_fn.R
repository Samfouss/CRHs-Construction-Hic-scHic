

create_bip_clust_graph_from_cell <- function(whole_cell_matrix, scHic_promoters_vec, cell_num, resolution = "6Mb"){
  
  par_default <- par(bty = 'n')
  row.names(whole_cell_matrix) = colnames(whole_cell_matrix)
  crhs_in_resolution = list()
  
  # Initilisation des des régions
  if(resolution == "1Mb"){
    resolution_part = c(1, 94, 187, 281, 375, 468, 562)
  }else if(resolution == "2Mb"){
    resolution_part = c(1, 187, 374, 562)
  }else if(resolution == "3Mb"){
    resolution_part = c(1, 281, 562)
  }else{
    resolution_part = c(1, 562)
  }
  
  for (res in seq_len(length(resolution_part)-1)) {
      
      cell_matrix = whole_cell_matrix[
        resolution_part[res]:resolution_part[res+1],
        resolution_part[res]:resolution_part[res+1],
        drop=FALSE
      ]
      
      promoters = scHic_promoters_ids[scHic_promoters_ids %in% resolution_part[res]:resolution_part[res+1]]
      
      promoters_names = row.names(whole_cell_matrix)[c(promoters)]
      promoters_ = which(row.names(cell_matrix) %in% promoters_names)
      
      
      if(length(promoters)>1){
        cell_matrix <- cell_matrix[
          promoters_, 
          -promoters_, 
          drop=FALSE
        ]
        
        
        net_bip <- graph_from_incidence_matrix(
          cell_matrix
        )
        
        if(!is_simple(net_bip)){
          net_bip <- simplify(net_bip, remove.multiple = FALSE, remove.loops = TRUE)
        }
        net_components_bip <- components(net_bip, mode = c("weak", "strong"))
        plot_main_cluster <- which(net_components_bip$csize>1)
        vert_ids <- V(net_bip)[net_components_bip$membership %in% plot_main_cluster]
        net_to_plot <- induced_subgraph(net_bip, vert_ids)
        V(net_to_plot)$color <- V(net_to_plot)$type
        V(net_to_plot)$color=gsub("FALSE","red",V(net_to_plot)$color)
        V(net_to_plot)$color=gsub("TRUE","lightblue",V(net_to_plot)$color)
        
        # A ce niveau, on prend chaque matrice d'incidence des différents CRHs qu'on ajoute à une liste 
        for (i in (1:net_components_bip$no)[net_components_bip$csize>1]){
          
          members <- net_components_bip$membership
          mat_bin = cell_matrix[
            rownames(cell_matrix) %in% names(members[members==i]),
            colnames(cell_matrix) %in% names(members[members==i]),
            drop=FALSE
          ]
          
          crhs_in_resolution[[length(crhs_in_resolution)+1]] = list(
            "name" = str_c("resolution.", resolution, ".", res),
            "mat_incidence" = mat_bin
          )
          
        }
        
        names(crhs_in_resolution) <- str_c("crh", 1:length(crhs_in_resolution))
        
      }
    }
  crhs_in_resolution
}


