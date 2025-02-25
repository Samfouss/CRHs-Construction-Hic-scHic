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

library(fpc)
library(igraph)
library(NbClust)
library(boot)

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

create_bip_clust_graph_from_cell <- function(whole_cell_matrix, scHic_promoters_vec, cell_num, resolution = "6Mb"){
  
  par_default <- par(bty = 'n')
  row.names(whole_cell_matrix) = colnames(whole_cell_matrix)
  crhs_in_resolution = list()
  
  # Initilisation des des régions
  if(resolution == "1Mb"){
    resolution_part = c(1, 94, 187, 281, 375, 468, 562)
  }else if(resolution == "2Mb"){
    resolution_part = c(1, 187, 374, 562)
  }else if(resolution == "3Mb"){
    resolution_part = c(1, 281, 562)
  }else{
    resolution_part = c(1, 562)
  }
  
  for (res in seq_len(length(resolution_part)-1)) {
    
    cell_matrix = whole_cell_matrix[
      resolution_part[res]:resolution_part[res+1],
      resolution_part[res]:resolution_part[res+1],
      drop=FALSE
    ]
    
    promoters = scHic_promoters_ids[scHic_promoters_ids %in% resolution_part[res]:resolution_part[res+1]]
    
    promoters_names = row.names(whole_cell_matrix)[c(promoters)]
    promoters_ = which(row.names(cell_matrix) %in% promoters_names)
    
    
    if(length(promoters)>1){
      cell_matrix <- cell_matrix[
        promoters_,
        -promoters_,
        drop=FALSE
      ]
      
      
      net_bip <- graph_from_incidence_matrix(
        cell_matrix
      )
      
      if(!is_simple(net_bip)){
        net_bip <- simplify(net_bip, remove.multiple = FALSE, remove.loops = TRUE)
      }
      net_components_bip <- components(net_bip, mode = c("weak", "strong"))
      plot_main_cluster <- which(net_components_bip$csize>1)
      vert_ids <- V(net_bip)[net_components_bip$membership %in% plot_main_cluster]
      net_to_plot <- induced_subgraph(net_bip, vert_ids)
      V(net_to_plot)$color <- V(net_to_plot)$type
      V(net_to_plot)$color=gsub("FALSE","red",V(net_to_plot)$color)
      V(net_to_plot)$color=gsub("TRUE","lightblue",V(net_to_plot)$color)
      
      # A ce niveau, on prend chaque matrice d'incidence des différents CRHs qu'on ajoute à une liste
      for (i in (1:net_components_bip$no)[net_components_bip$csize>1]){
        
        members <- net_components_bip$membership
        mat_bin = cell_matrix[
          rownames(cell_matrix) %in% names(members[members==i]),
          colnames(cell_matrix) %in% names(members[members==i]),
          drop=FALSE
        ]
        
        crhs_in_resolution[[length(crhs_in_resolution)+1]] = list(
          "name" = paste0("resolution.", resolution, ".", res),
          "mat_incidence" = mat_bin
        )
        
      }
      
      names(crhs_in_resolution) <- paste0("crh", 1:length(crhs_in_resolution))
      
    }
  }
  crhs_in_resolution
}
## Cette fonction permet de faire la comparaison entre les CRHs et donne les statistiques

# extend_matrix_with_names <- function(mat, mat_degeneration) {
#   r_1 = rownames(mat)
#   r_2 = rownames(mat_degeneration)
#   c_1 = colnames(mat)
#   c_2 = colnames(mat_degeneration)
#   # Calculate row and column names
#   rowmat <- union(r_1, r_2)
#   colmat <- union(c_1, c_2)
#
#   # Create resized matrices
#   mat_degeneration_redim <- matrix(-1, ncol = length(colmat), nrow = length(rowmat), dimnames = list(rowmat, colmat))
#   mat_redim <- mat_degeneration_redim
#
#   mat_redim[r_1, c_1] <- mat
#   mat_degeneration_redim[r_2, c_2] <- mat_degeneration
#
#   # Calculate sensitivity and specificity
#   sens <- sum(mat_degeneration_redim[mat_redim == mat_degeneration_redim] == 1) / sum(mat_degeneration_redim == 1)
#   spec <- sum(mat_degeneration_redim[mat_redim == mat_degeneration_redim] == 0) / sum(mat_degeneration_redim == 0)
#
#   return(c(sens, spec))
# }

compute_comparaison <- function(all_net_result, clu_chrs_result, make_degeneration = TRUE){
  
  ################ Création de la matrice servqnt à accuaillir les données de la comparaison ########
  
  ### Noms des lignes de la matrice : les CRHs puis leur block d'appartenance
  matLines = 0
  row_names = ""
  blocs = 2:length(all_net_result)
  for (l in blocs) {
    # print(length(all_net_result[[l]]$crhs))
    for (c in seq_len(length(all_net_result[[l]]$crhs))) {
      if(sum(all_net_result[[l]]$crhs[[c]]$mat_incidence) != -1){
        matLines = matLines +  1
        row_names = c(row_names, paste0("block_", l, "_crhs_", c))
      }
    }
  }
  row_names = row_names[row_names != ""]
  
  ### Noms des colonnes de la matrice : Les chrs puis leur cluster d'appartenance
  clus_num = c()
  for (cell in names(clu_chrs_result)) {
    n = as.numeric(unlist(strsplit(cell, "[_]"))[2])
    clus_num = c(clus_num, n)
  }
  matCol = 0
  col_names = ""
  for (clus in seq_len(length(clu_chrs_result))) {
    # print(length(clu_chrs_result[[clus]]$crhs))
    l = length(clu_chrs_result[[clus]])
    matCol = matCol + l
    col_names = c(col_names, paste0("cluster_", clus_num[clus], "_crhs_", seq_len(l)))
  }
  col_names = col_names[col_names != ""]
  
  ### Création de la matrice pour accueillir les résulats
  crhs_comparation_res = list(
    "sensibility_mat" = matrix(
      0,
      nrow = matLines,
      ncol = matCol,
      dimnames = list(row_names, col_names)
    ),
    "specificity_mat" = matrix(
      0,
      nrow = matLines,
      ncol = matCol,
      dimnames = list(row_names, col_names)
    )
  )
  
  col = 0
  
  # On parcours les cluster dans cette boucle
  for (clus in seq_along(clu_chrs_result)) {
    
    # On parcours les CRHs de chaque cluster dans cette boucle
    for (crh in seq_along(clu_chrs_result[[clus]])) {
      col = col + 1
      mat = clu_chrs_result[[clus]][[crh]]$mat_incidence
      # On reconstruit le réseau afin de recuperer après les arretes et noeuds à des fin de comparaison
      # net_bip_clus <- graph_from_incidence_matrix(mat)
      
      ln = 0
      # On boucle sur les blocks allant de 2 à 16
      for (bl in blocs) {
        block <- all_net_result[[bl]]
        crhs_count <- length(block$crhs)
        
        for (crh_ in seq_len(crhs_count)) {
          crh <- block$crhs[[crh_]]
          if (sum(crh$mat_incidence) != -1) {
            ln <- ln + 1
            mat_degeneration <- crh$mat_incidence
            if (make_degeneration) {
              mat_degeneration <- degenerationMatrix(mat_degeneration)
            }
            
            # redimenssion des deux matrices afin de calcluer les statistiques sur les intersections
            rowmat = union(rownames(mat), rownames(mat_degeneration))
            colmat = union(colnames(mat), colnames(mat_degeneration))
            
            mat_degeneration_redim <- matrix(
              -1,
              ncol=length(colmat),
              nrow=length(rowmat),
              dimnames=list(rowmat, colmat)
            )
            mat_redim <- mat_degeneration_redim
            
            indxA <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat), colnames(mat), FUN=paste)
            indxB <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat_degeneration), colnames(mat_degeneration), FUN=paste)
            mat_redim[indxA] <- mat
            mat_degeneration_redim[indxB] <- mat_degeneration
            
            crhs_comparation_res$sensibility_mat[ln, col] = sum(mat_degeneration_redim[mat_redim==mat_degeneration_redim]==1)/sum(mat_degeneration_redim==1)
            
            # mismatchs = length(rowmat[rowmat != rowmat_])*length(colmat[colmat != colmat_])
            crhs_comparation_res$specificity_mat[ln, col] = sum(mat_degeneration_redim[mat_redim==mat_degeneration_redim]==0)/sum(mat_degeneration_redim==0)
            }
          }
      }
    }
  }
  return(crhs_comparation_res)
}

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

set.seed(99999)
resolution = "3Mb"
load("rdata/umap.out_.rda")
load("rdata/cellUperDiagData.rda")
load("rdata/scHic_promoters_ids.rda")
load("rdata/all_net_result_complex_3Mb_.rda")

nb_cluster = 200:300
n <- 10000
umap_dt = umap.out$layout

all_net_result_complex_deg = all_net_result_complex_

###### Calcul du nombre de cluster optimal
compute_res = function(cells_clusters, cell_data){
  cluster_matrix_result = constr_mat_contacts(cells_clusters$cluster, cell_data)
  
  clu_chrs_result <- get_clusters_crhs(cluster_matrix_result, resolution = resolution)
  
  
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

# repositionner aléatoirement les lignes des donnnées afin que le choix des blocs ne soit pas lié à la position initiale
set.seed(123)
indices <- sample(n)

nb_echan = 20
k = 5

# Ici on procède à un shuffle des données
umap_dt = umap_dt[indices, ]
cellUperDiagData = cellUperDiagData[indices, ]


boot_fn = function(data, rows){
  # Ici on selectionne l'échantillon bootstrap
  umap_dt_ = umap_dt[rows, ]
  cellUperDiagData_ = cellUperDiagData[rows, ]
  
  # On procede à un shuffle pour la validation croisée
  indices_ <- sample(n)
  sen_vec = c(k)
  sep_vec = c(k)
  
  # Perform cross-validation
  for (i in 1:k) {
    # Calculate the start and end indices for the test set
    start <- ((i - 1) * (n / k)) + 1
    end <- i * (n / k)
    index <- indices_[start:end]
    # On forme nos données d'entrainement et de test ici
    test_data = umap_dt_[index, ]
    test_cell_data = cellUperDiagData_[index, ]
    
    train_data = umap_dt_[setdiff(indices_, index), ]
    train_cell_data = cellUperDiagData_[setdiff(indices_, index), ]
    
    rep_i <- NULL  # Initialize an empty variable to store the results
    
    for (clus in nb_cluster) {
      # Apply k-means clustering
      cells_clusters <- kmeans(train_data, clus, trace = FALSE)
      # Compute the results
      comp_res <- compute_res(cells_clusters, train_cell_data)
      # Append the results as a row
      rep_i <- rbind(rep_i, c(i, clus, comp_res[1], comp_res[2]))
    }
    
    # Determine optimal number of clusters based on sensitivity and specificity
    max_sen_indices <- which(rep_i[, 4] == max(rep_i[, 4]))
    if (length(max_sen_indices) > 1) {
      max_spec_indices <- max_sen_indices[which.max(rep_i[max_sen_indices, 3])]
      optimal_clus <- rep_i[max_spec_indices, 2]
    } else {
      optimal_clus <- rep_i[max_sen_indices, 2]
    }
    
    # Apply k-means on test data with optimal clusters and compute results
    cells_clusters <- kmeans(test_data, optimal_clus, trace = FALSE)
    comp_res <- compute_res(cells_clusters, test_cell_data)
    
    # Store sensitivity and specificity
    sep_vec[i] <- comp_res[1]
    sen_vec[i] <- comp_res[2]
    # res_nb_cluster = NbClust(data = train_data, distance = "euclidean", min.nc = 200, max.nc = 300, method = "kmeans")
    # save(res_nb_cluster, file = paste0("rdata/res_nb_cluster_", i,".rda"))
  }
  
  sep_vec_sum = summary(sep_vec)
  sen_vec_sum = summary(sen_vec)
  print("Fin pour un sous ensemble")
  return(c(sep_vec_sum[1], sep_vec_sum[2], sep_vec_sum[3], sep_vec_sum[4], sep_vec_sum[5], sep_vec_sum[6], sen_vec_sum[1], sen_vec_sum[2], sen_vec_sum[3], sen_vec_sum[4], sen_vec_sum[5], sen_vec_sum[6]))
}

r_boot = 1000
ncpus = 31

print("bootstrap debut")

boot_res <- boot(cellUperDiagData, boot_fn, R=r_boot, parallel = "multicore", ncpus = ncpus)

save(boot_res, file = "rdata/get_optimal_boot_res.rda")
print("End of bootstrap computation")