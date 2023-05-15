

source("scripts/crhs_creation_with_cells/cells_clustering_Kmeans.R")
source("scripts/crhs_creation_with_cells/juicerInputFileCreation_fn.R")

clusters_matrix <- list()
cell_dim = 562
clusters = unique(cells_clusters$cluster)

for (clus in seq_len(length(clusters))) {
  
  temp = cellUperDiagData[which(cells_clusters$cluster==clus), , drop = FALSE]
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
  
  juicerFile <- juicerInputCreation(mat_temp, ncells= 1)
  
  juicerFile <- juicerFile[juicerFile$score!= 0, ]
  
  write.table(
    juicerFile, 
    file = paste0("rdata/juicerInputFiles/cluster_", clus, "_juicer_input.txt"), 
    sep = "\t", 
    row.names = FALSE, 
    col.names = FALSE
  )
  
  mat_temp <- apply(mat_temp, 2, function(x) as.numeric(x>0))
  clusters_matrix[[length(clusters_matrix) + 1]] <- mat_temp
  
}

names(clusters_matrix) <- str_c("mat.inc.cluster_", clusters)


