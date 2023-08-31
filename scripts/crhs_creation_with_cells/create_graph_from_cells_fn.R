
create_bip_clust_graph_from_cell <- function(whole_cell_matrix, scHic_promoters_vec, cell_num, resolution = "6Mb"){
  
  par_default <- par(bty = 'n')
  row.names(whole_cell_matrix) = colnames(whole_cell_matrix)
  crhs_in_resolution = list()
  
  resolution_part = 562
  # Initilisation des des régions
  if(resolution == "1Mb"){
    resolution_part = c(1, 94, 187, 281, 375, 468, 562)
  }else if(resolution == "2Mb"){
    resolution_part = c(1, 187, 374, 562)
  }else if(resolution == "3Mb"){
    resolution_part = c(1, 281, 562)
  }
  
  if(length(resolution_part)==1){
    
    cell_matrix <- whole_cell_matrix[
      scHic_promoters_vec, 
      -scHic_promoters_vec, 
      drop=FALSE
    ]
    
    crhs = call_igraph_fn(cell_matrix, cell_num, resolution)
    
    crhs_in_resolution[[length(crhs_in_resolution)+1]] = crhs
    
  }else if(length(resolution_part)>1){
    for (res in seq_len(length(resolution_part)-1)) {
      
      cell_matrix = whole_cell_matrix[
        resolution_part[res]:resolution_part[res+1],
        resolution_part[res]:resolution_part[res+1],
        drop=FALSE
      ]
      
      promoters = scHic_promoters_ids[scHic_promoters_ids %in% resolution_part[res]:resolution_part[res+1]]
      promoters_names = row.names(whole_cell_matrix)[c(promoters)]
      promoters_ = which(row.names(cell_matrix) %in% promoters_names)
      
      # print(promoters_)

      if(length(promoters)>1){
        cell_matrix <- cell_matrix[
          promoters_, 
          -promoters_, 
          drop=FALSE
        ]
        # print(dim(cell_matrix))
        
        crhs = call_igraph_fn(cell_matrix, cell_num, resolution)
        # str(res)
        
        crhs_in_resolution[[length(crhs_in_resolution)+1]] = crhs
        
      }
    }
  }
  
  names(crhs_in_resolution) <- str_c("res", resolution, ".", seq_len(length(crhs_in_resolution)))
  crhs_in_resolution
}

