library("tidyverse")
library("FactoMineR")
library("factoextra")
library("ggpubr")

# construction de la matrice en retirant les CRHs non complexes
load("rdata/all_net_result.rda")
cells_matrix_ = list()

for (i in 2:length(all_net_result)) {
  
  crh_to_out = vector(mode = "numeric", length = length(all_net_result[[i]]$crhs))
  index = 0
  
  for (k in seq_len(length(all_net_result[[i]]$crhs))) {
    mat = all_net_result[[i]]$crhs[[k]]$mat_incidence
    
    if((nrow(mat) == 1 & ncol(mat) == 1) | (nrow(mat) == 2 & ncol(mat) == 1) | (nrow(mat) == 1 & ncol(mat) == 2)){
      crh_to_out[index] = k
      index = index + 1
    }
    
  }
  
  crh_to_out = crh_to_out[crh_to_out != 0]
  
  result <- all_net_result[[i]]$resume_fusion%>%
    str_split("-", simplify = TRUE)
  
  result <- gsub(",.*","",result)
  result <- gsub(".* ","",result)
  
  crhs_fus_len = nrow(result)
  
  compute_cells_matrix = matrix(
    0,
    nrow = 250,
    ncol = crhs_fus_len,
    dimnames = list(
      str_c("cell", sprintf("%03d", 1:250)),
      str_c("crhs", sprintf("%03d", 1:crhs_fus_len), "_block", sprintf("%02d", i))
    )
  )
  
  for(r in seq_len(crhs_fus_len)){
    
    # Pour une ligne donnée dans la metrice "result" des résultats, j'extrais les cellules conrespondant aux CRHs fusionnés
    extract_cells = str_sub(result[r, ][result[r, ] != ""], 1, 3)
    
    # Comme on ne compte qu'une seule fois les cellules meme si les deux chromosomes entrent dans la fusion, alors il nous suffit d'effacer les duplicate et de construire une matrice d'incidence à ajouter a la matrice finale "compute_cells_matrix"
    extract_cells = extract_cells[!duplicated(extract_cells)]
    
    compute_cells_matrix[
      str_c("cell", extract_cells),
      r
    ] <- 1
  }
  
  cells_matrix_[[length(cells_matrix_) + 1]] <- compute_cells_matrix[, -crh_to_out, drop = FALSE]
  
}


repetition_precess = function(all_net_result){
  
  for (i in 2:length(all_net_result)) {
    
    
    crh_to_out = vector(mode = "numeric", length = length(all_net_result[[i]]$crhs))
    index = 0
    
    for (k in seq_len(length(all_net_result[[i]]$crhs))) {
      if(length(all_net_result[[i]]$crhs[[k]])>1){
        mat = all_net_result[[i]]$crhs[[k]]$mat_incidence
        
        if((nrow(mat) == 1 & ncol(mat) == 1) | (nrow(mat) == 2 & ncol(mat) == 1) | (nrow(mat) == 1 & ncol(mat) == 2)){
          crh_to_out[index] = k
          index = index + 1
        }
      }
      
    }
    
    crh_to_out = crh_to_out[crh_to_out != 0]
    
    result <- all_net_result[[i]]$resume_fusion%>%
      str_split("-", simplify = TRUE)
    
    result <- gsub(",.*","",result)
    result <- gsub(".* ","",result)
    
    crhs_fus_len = nrow(result)
    
    compute_cells_matrix = matrix(
      0,
      nrow = 250,
      ncol = crhs_fus_len,
      dimnames = list(
        str_c("cell", sprintf("%03d", 1:250)),
        str_c("crhs", sprintf("%03d", 1:crhs_fus_len), "_block", sprintf("%02d", i))
      )
    )
    p
    for(r in seq_len(crhs_fus_len)){
      
      # Pour une ligne donnée dans la metrice "result" des résultats, j'extrais les cellules conrespondant aux CRHs fusionnés
      extract_cells = str_sub(result[r, ][result[r, ] != ""], 1, 3)
      
      # Comme on ne compte qu'une seule fois les cellules meme si les deux chromosomes entrent dans la fusion, alors il nous suffit d'effacer les duplicate et de construire une matrice d'incidence à ajouter a la matrice finale "compute_cells_matrix"
      extract_cells = extract_cells[!duplicated(extract_cells)]
      print(r)
      compute_cells_matrix[
        str_c("cell", extract_cells),
        r
      ] <- 1
    }
    
    cells_matrix_[[length(cells_matrix_) + 1]] <- compute_cells_matrix[, -crh_to_out, drop = FALSE]
    
  }
  
  cells_matrix_
}

repetition_precess(all_net_result_complex2_)

cells_matrix_ = repetition_precess(all_net_result)

cells_data_bind_ <- cells_matrix_[[1]]
for (b in 2:15) {
  cells_data_bind_ <- cbind(
    cells_data_bind_,
    cells_matrix_[[b]]
  )
}

save(cells_data_bind_, file = "rdata/all_rda_data/cells_data_bind_.rda")

dim(cells_data_bind_)
table(apply(cells_data_bind_, 2, sum))

# Ici on extrait les CRHs répétés parmi les CRHs assez complexes (cells_data_bind_)

# Répétition dans au moins 4 cellules
cells_data_4 <- cells_data_bind_[, apply(cells_data_bind_, 2, sum)>=4, drop = FALSE]
ncol(cells_data_4)
# 470
# Répétition dans au moins 5 cellules
cells_data_5 <- cells_data_bind_[, apply(cells_data_bind_, 2, sum)>=5, drop = FALSE]
ncol(cells_data_5)
# 174

cells_matrix_ = repetition_precess(all_net_result_complex2_)

cells_data_bind_ <- cells_matrix_[[1]]
for (b in 2:15) {
  cells_data_bind_ <- cbind(
    cells_data_bind_,
    cells_matrix_[[b]]
  )
}




