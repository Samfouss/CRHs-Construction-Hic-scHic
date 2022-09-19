# Chargement de packages
library("tidyverse")
library("FactoMineR")
library("factoextra")
library("ggpubr")

nb_class_max = 50
nb_class_min = 5

# construction de la matrice en retirant les CRHs non complexes
load("rdata/all_net_result.rda")
cells_matrix_ = list()

for (i in seq_len(length(all_net_result))) {
  
  crh_to_out = vector(mode = "numeric", length = length(all_net_result[[i]]$crhs))
  index = 0
  
  for (k in seq_len(length(all_net_result[[i]]$crhs))) {
    mat = all_net_result[[i]]$crhs[[k]]$mat_incidence
    
    if((nrow(mat) == 1 & ncol(mat) == 1) | (nrow(mat) == 2 & ncol(mat) == 1)){
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

cells_data_bind_ <- cells_matrix_[[1]]
for (b in 2:16) {
  cells_data_bind_ <- cbind(
    cells_data_bind_,
    cells_matrix_[[b]]
  )
}

# Clustering avec les données de la matrice construite en retirant les CRHs moins complexes

cells_data_ <- cells_data_bind_[, apply(cells_data_bind_, 2, sum)>2, drop = FALSE]
dim(cells_data_)

cells_data_[cells_data_ == 1] <- "in"
cells_data_[cells_data_ == 0] <- "Not in"

cells_data_ <- data.frame(cells_data_)

res.mca_ <- MCA(cells_data_, graph = FALSE)
eig.val_ <- get_eigenvalue(res.mca_)
fviz_mca_ind(res.mca_, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, ggtheme = theme_minimal())
which(res.mca_$eig[, 3, drop = FALSE]>=70)[1:5]

res.mca_ <- MCA(cells_data_, ncp = 249, graph = FALSE)
hcpc_cluster_ <- HCPC(res.mca_, kk=Inf, min = nb_class_min, max = nb_class_max, graph = FALSE)
fviz_cluster(hcpc_cluster_, repel = TRUE, geom = "point", main = "Classification des cellules après avoir retiré \n les CRHs moins complexes")



