# Créer deux matrices, une qu'on rempli par column et une par ligne

fill_symetric_matrix <- function(cells_clusters, cellUperDiagData){
  
  clusters_matrix <- list()
  cell_dim = 562
  clusters = unique(cells_clusters)
  for (clus in seq_len(length(clusters))) {
    
    temp = cellUperDiagData[which(cells_clusters==clus), , drop = FALSE]
    #temp <- apply(temp, 2, function(x) as.numeric(x>0))
    temp <- apply(temp, 2, sum)
    
    mat_row_name = paste0("BAC", sprintf("%03d", seq_len(cell_dim)))
    mat_temp <- matrix(
      0,
      nrow = cell_dim,
      ncol = cell_dim,
      dimnames = list(mat_row_name, mat_row_name)
    )
    
    mat_temp[upper.tri(mat_temp, diag = FALSE)] <- temp
    mat_temp <- mat_temp + t(mat_temp)
    
    mat_temp <- apply(mat_temp, 2, function(x) as.numeric(x>0))
    clusters_matrix[[length(clusters_matrix) + 1]] <- mat_temp
    
  }
  
  names(clusters_matrix) <- str_c("mat.inc.cluster_", clusters)
  clusters_matrix
}
