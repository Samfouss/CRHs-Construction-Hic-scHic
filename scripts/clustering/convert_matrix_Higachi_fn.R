library(vroom)

# To charge large dataset
# higachi_data <- fread("../Higashi/data/data.txt")
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
converting_matrix <- function(nb_cells, fenetre){
  
  # On charge une première matrice dont le nombre de lignes nous servira par la suite
  mat <- as.matrix(
    read.table(
      paste0("rdata/single_cell_hic_data/hic_mat_", sprintf("%03d", i), ".txt"), 
      quote="\"", 
      comment.char="", 
      stringsAsFactors = FALSE
    ) 
  )
  
  max_dim = nrow(mat)
  r = ceiling(max_dim/fenetre)
  # 
  if(max_dim<fenetre){
    
    cat("SVP vous devez donner une fênetre inférieur à la dimenssion de la matrice")
  
  }else{
    # 'cell_id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'count'
    cell_data_higachi <- matrix(
      NA,
      nrow = 2*r*nb_cells, 
      ncol = 6
    )
    
    for (i in seq_len(nb_cells)) {
      
      mat <- as.matrix(
        read.table(
          paste0("rdata/single_cell_hic_data/hic_mat_", sprintf("%03d", i), ".txt"), 
          quote="\"", 
          comment.char="", 
          stringsAsFactors = FALSE
        ) 
      )
      
      for (chr1 in (seq_len(r)-1)) {
        chr1_pos1 = chr1*fenetre + 1
        chr1_pos2 = (chr1+1)*fenetre
        if(chr1_pos2>max_dim){
          chr1_pos2 = max_dim
        }
        for (chr2 in (seq_len(r)-1)) {
          chr2_pos1 = chr2*fenetre + 1
          chr2_pos2 = (chr1+1)*fenetre
          if(chr2_pos2>max_dim){
            chr2_pos2 = max_dim
          }
          
          # cell_data_higachi[r*r*(i-1) + (r*chr1+(chr2 + 1)), 1] <- paste("cell_", i-1)
          # cell_data_higachi[r*r*(i-1) + (r*chr1+(chr2 + 1)), 2] <- i-1
          # cell_data_higachi[r*r*(i-1) + (r*chr1+(chr2 + 1)), 3] <- paste("chr", chr1+1)
          # cell_data_higachi[r*r*(i-1) + (r*chr1+(chr2 + 1)), 4] <- chr1_pos2
          # cell_data_higachi[r*r*(i-1) + (r*chr1+(chr2 + 1)), 5] <- paste("chr", chr1+1)
          # cell_data_higachi[r*r*(i-1) + (r*chr1+(chr2 + 1)), 6] <- 
          # cell_data_higachi[r*r*(i-1) + (r*chr1+(chr2 + 1)), 7] <- 
        }
      }
      
    } 
  }
  
}







