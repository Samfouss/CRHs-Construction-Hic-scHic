library(tidyverse)

load("rdata/all_net_result.rda")
result <- all_net_result$block1$resume_fusion%>%
  str_split("-", simplify = TRUE)

result <- gsub(",.*","",result)
result <- gsub(".* ","",result)

view(result)

crhs_fus_len = nrow(result)

cells_matrix = matrix(
  0,
  nrow = crhs_fus_len,
  ncol = 250,
  dimnames = list(
    str_c("crhs", sprintf("%03d", 1:crhs_fus_len)), 
    str_c("cell", sprintf("%03d", 1:250))
  )
)

for(r in seq_len(crhs_fus_len)){
  
  # Pour une ligne donnée dans la metrice "result" des résultats, j'extrais les cellules conrespondant aux CRHs fusionnés
  extract_cells = str_sub(result[r, ][result[r, ] != ""], 1, 3)
  
  # Comme on ne compte qu'une seule fois les cellules meme si les deux chromosomes entrent dans la fusion, alors il nous suffit d'effacer les duplicate et de construire une matrice d'incidence à ajouter a la matrice finale "cells_matrix"
  extract_cells = extract_cells[!duplicated(extract_cells)]
  
  cells_matrix[
    r,
    str_c("cell", extract_cells)
  ] <- 1
}

view(cells_matrix)

save(cells_matrix, file = "rdata/cells_matrix.rda")



