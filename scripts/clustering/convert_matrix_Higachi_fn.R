library(data.table)
library(tidyverse)

# To charge large dataset
higachi_data <- fread("../Higashi/data/data.txt")
# nb_cells = 10

# Cette onction est créée pour convertir le fichier des données scHic simulées en fichier accepté par le programme HIgachi
# Le fichier de donné doit être structuré de la sorte 'cell_id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'count'
# fenetre : c'est le nombre de billes que l'on veut considérer pour sommer les contacts de la matrice
# nb_cells : le nombre de cellules à lire

##### Forat du fichier a obtenir pour higachi
#     cell_name             cell_id chrom1      pos1 chrom2      pos2 count
# 1: ML1_GAGGAGCA_CGATGACA       0   chr2  88637036   chr2  89109729     1
# 2: ML1_GCTACGGT_AGTCGTAT       1   chr4 126924974   chr4 127334997     1
# 3: ML1_AGGTGCGA_ATACATGT       2  chr15  42130039  chr15  42677209     1
# 4: ML1_GCCTCGAA_GAGTACGT       3   chr1 209393848   chr1 232468122     1
# 5: ML1_GCTCGCTA_CTAGTGAA       4   chr5   5565796   chr5 149955845     1
# 6: ML1_CAGGCTTG_GATATAAC       5   chr1 181632778   chr1 181731968     1

# Une bille fait 2667 bases
converting_matrix <- function(nb_cells, nb_mat_line){
  # Nombre de bases que représente une bille dans les polymeres
  bins_equiv = 2667

  cell_dim = nb_mat_line*(nb_mat_line-1)/2
  nb_rows = cell_dim*nb_cells

  # 'cell_id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'count'
  cell_data_higachi <- matrix(
    NA,
    nrow = nb_rows, 
    ncol = 6
  )
  
  #cell_data_higachi = as.data.frame(cell_data_higachi)
  
  for (i in seq_len(nb_cells)) {
    
    mat <- as.matrix(
      read.table(
        paste0("rdata/single_cell_hic_data/hic_mat_", sprintf("%03d", i), ".txt"), 
        quote="\"", 
        comment.char="", 
        stringsAsFactors = FALSE
      ) 
    )
    
    block_cell = cell_dim*(i-1)
    # Initialisation de l'index de la ligne
    line = 0
    
    for (chr1 in seq_len(nb_mat_line)) {
      # pos1 = ceiling(((chr1-1)*bins_equiv + chr1*bins_equiv)/2)
      pos1 = ceiling(bins_equiv*(chr1-1) + bins_equiv/2)
      
      for (chr2 in seq_len(nb_mat_line)) {
        # Ici, il nous faut juste prendre la partie supérieure (ou inférieure), la diagonale exlue de la matrice des contacts
        # Donc pour une ligne donnée, on prend les céllules dont l'index de la colonne est superieur à celui de la ligne en question
        if(chr2 > chr1){
          pos2 = ceiling(bins_equiv*(chr2-1) + bins_equiv/2)
          line = line + 1

          cell_data_higachi[block_cell+line, 1] <- i-1
          cell_data_higachi[block_cell+line, 2] <- "chr1"
          cell_data_higachi[block_cell+line, 3] <- pos1
          cell_data_higachi[block_cell+line, 4] <- "chr2"
          cell_data_higachi[block_cell+line, 5] <- pos2
          cell_data_higachi[block_cell+line, 6] <- mat[chr1, chr2]
        }
      }
    }
    
  }

  colnames(cell_data_higachi) <- c('cell_id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'count')
  cell_data_higachi
  
}


cell_data_higachi <- converting_matrix(100, 562)

# Save promoters ID
save(cell_data_higachi, file = "rdata/cell_data_higachi.rda")

# write.table(cell_data_higachi, file = "rdata/data.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

