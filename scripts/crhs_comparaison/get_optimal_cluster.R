# Function to check and install a package if not installed
install_if_not_installed <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

# Example: Checking and installing the 'ggplot2' package
install_if_not_installed("tidyverse")
install_if_not_installed("igraph")
install_if_not_installed("fpc")
install_if_not_installed("NbClust")
install_if_not_installed("doParallel")
install_if_not_installed("foreach")

library(fpc)
library(igraph)
library(NbClust)
library(boot)
library(doParallel)
library(foreach)

set.seed(99999)
resolution = "3Mb"
load("rdata/umap.out_.rda")
load("rdata/cellUperDiagData.rda")
load("rdata/scHic_promoters_ids.rda")
load("rdata/all_net_result_complex_3Mb_.rda")

nb_cluster = 200:300
n <- 10000
umap_dt = umap.out$layout
k = 5

res_matrix = matrix(NA, nrow =2 , ncol = 4)

# repositionner aléatoirement les lignes des donnnées afin que le choix des blocs ne soit pas lié à la position initiale
indices <- sample(n)

nb_prom = c()
nb_enhan = c()

# Ici on procède à un shuffle des données
umap_dt = umap_dt[indices, ]
cellUperDiagData = cellUperDiagData[indices, ]

# Set up parallel backend to use available cores
num_cores <- detectCores() - 1  # Leave one core free for other tasks
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# On procede à un shuffle pour la validation croisée
indices_ <- sample(n)
sen_vec = c()
sep_vec = c()

##### Fonction permettant de construire les matrices d'incidence
constr_mat_contacts <- function(cells_clusters, cellUperDiagData, get_incidence_mat=TRUE){
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
    if (get_incidence_mat){
      mat_temp <- apply(mat_temp, 2, function(x) as.numeric(x>0))
    }
    clusters_matrix[[length(clusters_matrix) + 1]] <- mat_temp
    
  }
  if (get_incidence_mat){
    names(clusters_matrix) <- paste0("mat.inc.cluster_", clusters)
  }else{
    names(clusters_matrix) <- paste0("mat.cluster_", clusters)
  }
  clusters_matrix
}

source("create_graph_from_cells_fn.R")
## Cette fonction permet de faire la comparaison entre les CRHs et donne les statistiques

source("compute_comparaison.R")

degenerationMatrix <- function(matrixToDegenerate){
  
  # matrixToDegenerate = matrixToDegenerate
  
  # Extraction des noms des lignes et colonnes de la matrice à degenerer
  row_names = row.names(matrixToDegenerate)
  col_names = colnames(matrixToDegenerate)
  # Recupération des numeros de blocks (inutile)
  # bins_color = unique(str_sub(colnames(matrixToDegenerate), 1, 3))
  # Récupération du numero de la bille. Ces numero seront utilisés pour identifier la position et le bac de chaque bille
  row_names_numb = substr(row_names, start = 4, stop = 7)
  col_names_numb = substr(col_names, start = 4, stop = 7)
  
  # A partir de la matrice de transition, on identifie le bac d'appartenance de chaque bille en fonction de sa position
  row_names_ = unique(bacs_matrix[which(bacs_matrix[, 1]%in%row_names_numb), 2])
  col_names_ = unique(bacs_matrix[which(bacs_matrix[, 1]%in%col_names_numb), 2])
  
  # Initilisation de la matrice de réduction
  mat_reduct <- matrix(
    0,
    nrow = length(row_names_),
    ncol = length(col_names_),
    dimnames = list(
      row_names_,
      col_names_
    )
  )
  
  for (bac_i in row_names_) {
    for (bac_j in col_names_) {
      # On recupere ici les indices en ligne de toutes les billes faisant partie du bac i
      i = which(substr(row_names, start = 4, stop = 7) %in% bacs_matrix[bacs_matrix[, 2]==bac_i, 1])
      # On recupere ici les indices en ligne de toutes les billes faisant partie du bac j
      j = which(substr(col_names, start = 4, stop = 7) %in% bacs_matrix[bacs_matrix[, 2]==bac_j, 1])
      # On sum les contacts entre les billes recuperées en ligne et en colonne
      mat_reduct[bac_i, bac_j] = 0
      if(sum(matrixToDegenerate[i, j]) > 0){
        mat_reduct[bac_i, bac_j] = 1
      }
    }
  }
  
  # On retourne la matrice construite
  mat_reduct
  
}
################## Création de la matrice de transition permettant de passer des billes au bacs ######
bins_number = 1:2250
bacs_number = c(0, which(bins_number%%4==0))
bacs_matrix = matrix(
  "",
  nrow = length(bins_number),
  ncol = 2
)

for (line in seq_len(length(bins_number))) {
  b = which(bacs_number == (line-line%%4))
  if(line%%4==0){
    b = which(bacs_number == (line-line%%4)) - 1
  }
  bacs_matrix[line, 1] = sprintf("%04d", line)
  bacs_matrix[line, 2] = paste0("BAC", sprintf("%03d", b))
}

rm("bins_number", "bacs_number", "line", "b")

get_clusters_crhs <- function(clusters_matrix, resolution = "6Mb"){
  clu_chrs_result = list()
  for (clus in seq_len(length(clusters_matrix))) {
    net = create_bip_clust_graph_from_cell(
      clusters_matrix[[clus]],
      scHic_promoters_ids,
      clus,
      resolution
    )
    clu_chrs_result[[length(clu_chrs_result)+1]] <- net
  }
  
  clu_chrs_result
}
################################################################################

all_net_result_complex_deg = all_net_result_complex_
for (bl in 2:16) {
  for (i in seq_len(length(all_net_result_complex_[[bl]]$crhs))) {
    mat = all_net_result_complex_[[bl]]$crhs[[i]]$mat_incidence
    if(sum(mat)!=-1){
      all_net_result_complex_deg[[bl]]$crhs[[i]][[2]]= degenerationMatrix(mat)
    }
  }
}

###### Calcul du nombre de cluster optimal
compute_res = function(cells_clusters, cell_data){
  cluster_matrix_result = constr_mat_contacts(cells_clusters$cluster, cell_data)
  
  clu_chrs_result <- get_clusters_crhs(cluster_matrix_result, resolution = resolution)
  
  for (i in seq_len(length(clu_chrs_result))) {
    for (j in seq_len(length(clu_chrs_result[[i]]))) {
      nb_prom = c(nb_prom, dim(clu_chrs_result[[i]][[j]]$mat_incidence)[1])
      nb_enhan = c(nb_enhan, dim(clu_chrs_result[[i]][[j]]$mat_incidence)[2])
    }
  }
  
  crhs_comparaison_res = compute_comparaison(all_net_result_complex_deg, clu_chrs_result, make_degeneration = FALSE)
  
  nb_crhs_in_clus = dim(crhs_comparaison_res$specificity_mat)[2]
  nb_crh_in_struc = dim(crhs_comparaison_res$specificity_mat)[1]
  
  highest_spec <- matrix(
    "",
    nrow = nb_crh_in_struc,
    ncol = 3,
  )
  highest_sen <- matrix(
    "",
    nrow = nb_crh_in_struc,
    ncol = 3,
  )
  
  sen_mat = crhs_comparaison_res$sensibility_mat
  spec_mat = crhs_comparaison_res$specificity_mat
  
  l1 = 1
  l2 = 1
  for (res in seq_len(nb_crh_in_struc)){
    if(!all(is.na(spec_mat[res, ]))){
      highest_spec[l1, 1] <- names(spec_mat[, 1][res])
      highest_spec[l1, 2] <- names(which.max(spec_mat[res, ]))
      highest_spec[l1, 3] <- max(spec_mat[res, ][!is.na(spec_mat[res, ])])
      l1 = l1 + 1
    }
    if(!all(is.na(sen_mat[res, ]))){
      highest_sen[l2, 1] <- names(sen_mat[, 1][res])
      highest_sen[l2, 2] <- names(which.max(sen_mat[res, ]))
      highest_sen[l2, 3] <- max(sen_mat[res, ][!is.na(sen_mat[res, ])])
      l2 = l2 + 1
    }
  }
  
  highest_spec = highest_spec[highest_spec[, 1] != "", ]
  highest_spec_ = as.numeric(highest_spec[, 3])
  
  highest_sen = highest_sen[highest_sen[, 1] != "", ]
  highest_sen_ = as.numeric(highest_sen[, 3])
  
  return(c(mean(highest_spec_), mean(highest_sen_)))
}

# Perform cross-validation
for (i in 1:k) {
  print(paste0(i, "-", format(Sys.time(), "%H:%M:%S")))
  # Calculate the start and end indices for the test set
  start <- ((i - 1) * (n / k)) + 1
  end <- i * (n / k)
  index <- indices_[start:end]
  
  # Split data into training and testing sets
  test_data <- umap_dt[index, ]
  test_cell_data <- cellUperDiagData[index, ]
  
  train_data <- umap_dt[setdiff(indices_, index), ]
  train_cell_data <- cellUperDiagData[setdiff(indices_, index), ]
  
  # Run clustering in parallel over nb_cluster
  print(format(Sys.time(), "%H:%M:%S"))
  rep_i <- foreach(clus = nb_cluster, .combine = rbind) %dopar% {
    cells_clusters <- kmeans(train_data, clus, trace = FALSE)
    comp_res <- compute_res(cells_clusters, train_cell_data)
    c(i, clus, comp_res[1], comp_res[2])
  }
  print(format(Sys.time(), "%H:%M:%S"))
  
  # Find the optimal number of clusters based on sensitivity and specificity
  max_sen <- which(rep_i[, 4] == max(rep_i[, 4]))
  
  if (length(max_sen) > 1) {
    max_spec <- which(rep_i[max_sen, 3] == max(rep_i[max_sen, 3]))
    if (length(max_spec) > 1) {
      optimal_clus <- min(rep_i[max_spec, 2])
    } else {
      optimal_clus <- rep_i[max_spec, 2]
    }
  } else {
    optimal_clus <- rep_i[max_sen, 2]
  }
  
  # Perform k-means on the test set with the optimal number of clusters
  cells_clusters <- kmeans(test_data, optimal_clus, trace = FALSE)
  comp_res <- compute_res(cells_clusters, test_cell_data)
  
  # Collect the results for this fold
  res_matrix <- rbind(res_matrix, c(i, optimal_clus, comp_res[1], comp_res[2]))
  sep_vec <- c(sep_vec, comp_res[1])
  sen_vec <- c(sen_vec, comp_res[2])
  
  # Optional: Save results per fold if needed
  # res_nb_cluster = NbClust(data = train_data, distance = "euclidean", min.nc = 200, max.nc = 300, method = "kmeans")
  # save(res_nb_cluster, file = paste0("rdata/res_nb_cluster_", i,".rda"))
  print(paste0("Fin ", i, format(Sys.time(), "%H:%M:%S")))
}

# Stop the parallel backend
stopCluster(cl)

save(sep_vec, file = "rdata/sep_vec2.rda")
save(sen_vec, file = "rdata/sen_vec2.rda")

save(nb_prom, file = "rdata/nb_prom2.rda")
save(nb_enhan, file = "rdata/nb_enhan2.rda")

print("End of bootstrap computation")

save(res_matrix, file = "rdata/get_optimal_res_matrix2.rda")
print("End")