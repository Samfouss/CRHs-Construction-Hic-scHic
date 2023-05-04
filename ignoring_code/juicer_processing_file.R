library(data.table)
library(tidyverse)


start_ = 109000000
bins = 2667
bins_per_bacs = 4
mat_dim = 562
df_nrow = mat_dim*(mat_dim-1)/2

# Une bille fait 2667 bases
converting_into_short_format <- function(nb_cells){
  
  for (i in seq_len(nb_cells)) {
    
    # mat <- as.matrix(
    #   read.table(
    #     paste0("rdata/single_cell_hic_data/hic_mat_", sprintf("%03d", i), ".txt"), 
    #     quote="\"", 
    #     comment.char="", 
    #     stringsAsFactors = FALSE
    #   ) 
    # )
    
    # <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2>
    
    cell_data <- data.frame(
      str1 = rep(0, df_nrow),
      chr1 = rep("chr11", df_nrow),
      pos1 = c(NA),
      frag1 = rep(0, df_nrow),
      str2 = rep(0, df_nrow),
      chr2 = rep("chr11", df_nrow),
      pos2 = c(NA),
      frag2 = rep(1, df_nrow)
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
          
        }
      }
    }
    
  }
  
  write.table(
    cell_data, 
    file = paste0("rdata/single_cell_hic_data/juicer_input_", sprintf("%03d", i), ".txt"), 
    sep = "\t", 
    row.names = FALSE, 
    col.names = FALSE
  )
  
}


converting_into_short_format(250)


