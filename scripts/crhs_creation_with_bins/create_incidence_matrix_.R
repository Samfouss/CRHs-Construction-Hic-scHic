library(tidyverse)
load("rdata/all_rda_data/promoters_ids.rda")

ncells = 250
min_dist = 3


for (cell in seq_len(ncells)) {
  
  data1 = read.table(paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/paire_1.txt"))%>%
    mutate(
      ID = paste0("B", sprintf("%02d", V4), sprintf("%04d", 1:n()))
    )
  
  data2 = read.table(paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/paire_2.txt"))%>%
    mutate(
      ID = paste0("B", sprintf("%02d", V4), sprintf("%04d", 1:n()))
    )
  
  names1 = data1[, "ID"]
  names2 = data2[, "ID"]

  compute_dist1 <- matrix(
    0,
    nrow = length(names1),
    ncol = length(names1),
    dimnames = list(names1, names1)
  )
  
  compute_dist2 <- matrix(
    0,
    nrow = length(names2),
    ncol = length(names2),
    dimnames = list(names2, names2)
  )
  
  for (bl in 2:16) {
    
    mat_row_name1 = data1[data1$V4==bl, "ID"]
    dist_data1 = data1[data1$V4==bl, 1:3]

    compute_dist1_ <- as.matrix(
      dist(
        x = dist_data1, 
        method = "euclidean", 
        diag = TRUE
      )
    )
    
    compute_dist_bin1 <- sapply(
      compute_dist1_, 
      FUN = function(item) ifelse(item<=min_dist, 1, 0)
    )
    
    compute_dist_bin1 <- matrix(
      compute_dist_bin1, 
      nrow = length(mat_row_name1), 
      ncol = length(mat_row_name1), 
      byrow = TRUE,
      dimnames = list(mat_row_name1, mat_row_name1)
    )
    
    mat_row_name2 = data2[data2$V4==bl, "ID"]
    dist_data2 = data2[data2$V4==bl, 1:3]
    
    compute_dist2_ <- as.matrix(
      dist(
        x = dist_data2, 
        method = "euclidean", 
        diag = TRUE
      )
    )
    
    compute_dist_bin2 <- sapply(
      compute_dist2_, 
      FUN = function(item) ifelse(item<=min_dist, 1, 0)
    )
    
    compute_dist_bin2 <- matrix(
      compute_dist_bin2, 
      nrow = length(mat_row_name2), 
      ncol = length(mat_row_name2), 
      byrow = TRUE,
      dimnames = list(mat_row_name2, mat_row_name2)
    )
    
    # print(cell)
    # print(bl)
    # print(mat_row_name1)
    
    compute_dist1[mat_row_name1, mat_row_name1] = compute_dist_bin1[mat_row_name1, mat_row_name1]
    compute_dist2[mat_row_name2, mat_row_name2] = compute_dist_bin2[mat_row_name2, mat_row_name2]
    
  }
  
  compute_dist2 = compute_dist2[promoters_ids, -promoters_ids]
  compute_dist1 = compute_dist1[promoters_ids, -promoters_ids]
  
  write.table(compute_dist1, file=paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/incidence_matrix_1.txt"), row.names=TRUE, col.names=TRUE)
  write.table(compute_dist2, file=paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/incidence_matrix_2.txt"), row.names=TRUE, col.names=TRUE)
  
  
}




for (cell in seq_len(ncells)) {
  
  mat1 = read.table(paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/incidence_matrix_1.txt"))
  
  mat2 = read.table(paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/incidence_matrix_2.txt"))
  
  rowsn = union(row.names(mat1), row.names(mat2))
  coln = union(colnames(mat1), colnames(mat2))
  
  compute_dist <- matrix(
    0,
    nrow = length(rowsn),
    ncol = length(coln),
    dimnames = list(rowsn, coln)
  )
  
  indxA <- outer(rowsn, coln, FUN=paste) %in% outer(rownames(mat1), colnames(mat1), FUN=paste)
  indxB <- outer(rowsn, coln, FUN=paste) %in% outer(rownames(mat2), colnames(mat2), FUN=paste)
  # On récupère ici les données dans les matrices redimenssionnées. Avec les noms des lignes et colonnes (indexes), on est sur de recuprerer les bonnes données
  compute_dist[indxA] <- mat1
  compute_dist[indxB] <- mat2
  
  write.table(compute_dist, file=paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/incidence_matrix.txt"), row.names=TRUE, col.names=TRUE)
  
  print(cell)
}





