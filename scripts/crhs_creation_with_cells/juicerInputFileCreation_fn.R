
################ Fonction permettant de produire des table de données adaptées à Juicer ##############

juicerInputCreation <- function(matrix_to_process, start_ = 109000000, bins = 2667, bins_per_bacs = 4, mat_dim = 562, ncells = 1){
  
  df_nrow = mat_dim*(mat_dim-1)/2
  
  cell_data <- data.frame(
    str1 = rep(0, df_nrow),
    chr1 = rep(11, df_nrow),
    pos1 = c(NA),
    frag1 = rep(0, df_nrow),
    str2 = rep(0, df_nrow),
    chr2 = rep(11, df_nrow),
    pos2 = c(NA),
    frag2 = rep(1, df_nrow),
    score = rep(0, df_nrow)
  )
  
  for (bac_i in seq_len(mat_dim)) {
    
    bac_start = start_ + (bac_i-1)*bins*bins_per_bacs
    bac_end = start_ + bac_i*bins*bins_per_bacs
    pos1 = ceiling((bac_start + bac_end)/2)
    
    for (bac_j in seq_len(mat_dim)) {
      
      if(bac_j > bac_i){
        
        bac_start = start_ + (bac_j-1)*bins*bins_per_bacs
        bac_end = start_ + bac_j*bins*bins_per_bacs
        pos2 = ceiling((bac_start + bac_end)/2)
        
        first_NA_line = which(is.na(cell_data$pos1))[1]
        cell_data[first_NA_line, "pos1"] <- pos1
        cell_data[first_NA_line, "pos2"] <- pos2
        cell_data[first_NA_line, "score"] <- round(matrix_to_process[bac_i, bac_j]/ncells, 4)
        
      }
    }
  }
  
  cell_data
  
}



####### Construction des matrices de contacts par cluster #############

constr_mat_contacts <- function(cells_clusters, cellUperDiagData){
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
    # cells_clusters$size[clus]
    juicerFile <- juicerInputCreation(mat_temp, ncells= 1)
    juicerFile$score = round(juicerFile$score)
    
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
  clusters_matrix
}

