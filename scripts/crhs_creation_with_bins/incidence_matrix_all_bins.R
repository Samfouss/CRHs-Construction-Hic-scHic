

# Chargement des données
source("./scripts/load_save_data.R")

# Recupération du nombre de structure
nb_replicas = length(unique(all_paired_structure$paire))

nb_replicas = 1
# La boucle ira de 1 à 500
#nb_replicas = 50
for (r in 1:nb_replicas) {
  
  # Construction de la matrice d'incidence avec la première structure de la cellule
  str = 1
  data <- all_paired_structure%>%
    filter(paire == str_c(r, sprintf("%03d", str)))%>%
    select(-c(ends_with("_c"), "paire"))%>%
    mutate(
      ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n()))
    )
  
  three_dim_data <- data[, 1:3]
  mat_row_name = data[, 5]
  compute_dist <- as.matrix(
    dist(
      x = three_dim_data, 
      method = "euclidean", 
      diag = TRUE
    )
  )
  
  compute_dist_bin1 <- sapply(
    compute_dist, 
    FUN = function(item) ifelse(item<=3, 1, 0)
  )
  
  compute_dist_bin1 <- matrix(
    compute_dist_bin1, 
    nrow = nrow(compute_dist), 
    ncol = ncol(compute_dist), 
    byrow = TRUE,
    dimnames = list(mat_row_name, mat_row_name)
  )
  
  diag(compute_dist_bin1) <- 0
  
  # compute_dist_bin1 est la première matrice d'incidence. Il faut mettre à présent à 0 toutes les interaction entre les billes de différentes couleures de la matrice. Il nous suffit donc de verifier si les billes en ligne et en colonne sont du meme bloc
  for (l in mat_row_name) {
    for (c in mat_row_name) {
      if (str_sub(l, 1, 3) != str_sub(c, 1, 3)){
        compute_dist_bin1[l, c] = 0
      }
    }
  }
}

# Nombre de contact trouvé dans la matrice d'incidence
sum(compute_dist_bin1)/2
# La matrice est bien symétrique
isSymmetric(compute_dist_bin1)

# Le nombre de contacts dans la matrice simulée
matrix_sum = rep(0, 10)
for (i in 1:10) {
  
  matrix_sum[i] = sum(
    read.table(paste0("Sc_Cell_In-silico_Hi-C/hic_mat_", sprintf("%03d", i), ".txt"), quote="\"", comment.char="", stringsAsFactors = FALSE)
  )
  
}

# La proportion du nombre de contact recupéré
2706/sum(compute_dist_bin1)








