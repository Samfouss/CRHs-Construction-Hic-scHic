library(tidyverse)
load("rdata/all_rda_data/promoters_ids.rda")

ncells = 250
min_dist = 3

### Production of incidence matrix for each chromosome
for (cell in seq_len(ncells)) {
  
  # read data for chromosome 1
  data1 = read.table(paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/paire_1.txt"))%>%
    mutate(
      ID = paste0("B", sprintf("%02d", V4), sprintf("%04d", 1:n()))
    )
  # read data for  chromosome 2
  data2 = read.table(paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/paire_2.txt"))%>%
    mutate(
      ID = paste0("B", sprintf("%02d", V4), sprintf("%04d", 1:n()))
    )
  
  # Get just IDs of each data
  names1 = data1[, "ID"]
  names2 = data2[, "ID"]

  # Create matrix which gonna contain the contacts for data 1. So for all blocks and beads
  compute_dist1 <- matrix(
    0,
    nrow = length(names1),
    ncol = length(names1),
    dimnames = list(names1, names1)
  )
  
  # Create matrix which gonna contain the contacts for data 2. So for all blocks and beads
  compute_dist2 <- matrix(
    0,
    nrow = length(names2),
    ncol = length(names2),
    dimnames = list(names2, names2)
  )
  
  # Loop for each block in the data
  for (bl in 2:16) {
    
    ########### For the first matrix
    # Get ID of beads in the current block
    mat_row_name1 = data1[data1$V4==bl, "ID"]
    # Get the tree dimensional data of beads in the current block
    dist_data1 = data1[data1$V4==bl, 1:3]
    # Compute 2D distance matrix with bead in the current block
    compute_dist1_ <- as.matrix(
      dist(
        x = dist_data1, 
        method = "euclidean", 
        diag = TRUE
      )
    )
    # Now get just the one which are closer (distance less than 3). 
    compute_dist_bin1 <- sapply(
      compute_dist1_, 
      FUN = function(item) ifelse(item<=min_dist, 1, 0)
    )
    # Beads which are closer get 1 in this matrix and 0 if not
    compute_dist_bin1 <- matrix(
      compute_dist_bin1, 
      nrow = length(mat_row_name1), 
      ncol = length(mat_row_name1), 
      byrow = TRUE,
      dimnames = list(mat_row_name1, mat_row_name1)
    )
    ################ Repaet the same process (like previous) for the second chromosome 
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
    
    # Once we get the incidence matrix for each chromosome, we gonna put data in the big matrix, the one which gonna contain all data
    compute_dist1[mat_row_name1, mat_row_name1] = compute_dist_bin1[mat_row_name1, mat_row_name1]
    compute_dist2[mat_row_name2, mat_row_name2] = compute_dist_bin2[mat_row_name2, mat_row_name2]
    
  }
  # One we get all data, here we remove the promoters in column and just keep them in rows. We do the same thing with the two chromosome
  compute_dist2 = compute_dist2[promoters_ids, -promoters_ids]
  compute_dist1 = compute_dist1[promoters_ids, -promoters_ids]
  
  # Save data
  write.table(compute_dist1, file=paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/incidence_matrix_1.txt"), row.names=TRUE, col.names=TRUE)
  write.table(compute_dist2, file=paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/incidence_matrix_2.txt"), row.names=TRUE, col.names=TRUE)
  
}


### Production of incidence matrix for the cell
for (cell in seq_len(ncells)) {
  # load matrix data
  mat1 = as.matrix(read.table(paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/incidence_matrix_1.txt")))
  mat2 = as.matrix(read.table(paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/incidence_matrix_2.txt")))
  # Get row and column rows. Here we get the unique beads IDs from each matrix by doing a union.
  rowsn = union(row.names(mat1), row.names(mat2))
  coln = union(colnames(mat1), colnames(mat2))
  # Create matrix which gonna contain the union data for the first data
  compute_dist_ <- matrix(
    0,
    nrow = length(rowsn),
    ncol = length(coln),
    dimnames = list(rowsn, coln)
  )
  # Create matrix which gonna contain the union data for the second data
  compute_dist__ = compute_dist_
  # Get index for the each matrix to get data
  indxA <- outer(rowsn, coln, FUN=paste) %in% outer(rownames(mat1), colnames(mat1), FUN=paste)
  indxB <- outer(rowsn, coln, FUN=paste) %in% outer(rownames(mat2), colnames(mat2), FUN=paste)
  
  # get data for each chromosom in the union matrix
  compute_dist_[indxA] = mat1
  compute_dist__[indxB] = mat2
  # Here we gonna sum the matrix. In that way we gonna get all data
  compute_dist = compute_dist_ + compute_dist__
  compute_dist[compute_dist>0] = 1
  # Save data
  write.table(compute_dist, file=paste0("rdata/cell_folder/cell_", sprintf("%03d", cell), "/incidence_matrix.txt"), row.names=TRUE, col.names=TRUE)
  
  # print(cell)
}





