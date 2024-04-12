# Chargement des librairies
library(fpc)
library(tidyverse)
library(igraph)
library(boot)
library(umap)
library(ggplot2)

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
    names(clusters_matrix) <- str_c("mat.inc.cluster_", clusters)
  }else{
    names(clusters_matrix) <- str_c("mat.cluster_", clusters)
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
          "name" = str_c("resolution.", resolution, ".", res),
          "mat_incidence" = mat_bin
        )
        
      }
      
      names(crhs_in_resolution) <- str_c("crh", 1:length(crhs_in_resolution))
      
    }
  }
  crhs_in_resolution
}

## Cette fonction permet de faire la comparaison entre les CRHs et donne les statistiques

extend_matrix_with_names <- function(mat, mat_degeneration) {
  r_1 = rownames(mat)
  r_2 = rownames(mat_degeneration)
  c_1 = colnames(mat)
  c_2 = colnames(mat_degeneration)
  # Calculate row and column names
  rowmat <- union(r_1, r_2)
  colmat <- union(c_1, c_2)
  
  # Create resized matrices
  mat_degeneration_redim <- matrix(-1, ncol = length(colmat), nrow = length(rowmat), dimnames = list(rowmat, colmat))
  mat_redim <- mat_degeneration_redim
  
  mat_redim[r_1, c_1] <- mat
  mat_degeneration_redim[r_2, c_2] <- mat_degeneration
  
  # Calculate sensitivity and specificity
  sens <- sum(mat_degeneration_redim[mat_redim == mat_degeneration_redim] == 1) / sum(mat_degeneration_redim == 1)
  spec <- sum(mat_degeneration_redim[mat_redim == mat_degeneration_redim] == 0) / sum(mat_degeneration_redim == 0)
  
  return(c(sens, spec))
}

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
        row_names = c(row_names, str_c("block_", l, "_crhs_", c))
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
    col_names = c(col_names, str_c("cluster_", clus_num[clus], "_crhs_", seq_len(l)))
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
            
            res = extend_matrix_with_names(mat, mat_degeneration)
            
            # Assign to result matrices
            crhs_comparation_res$sensibility_mat[ln, col] <- res[1]
            crhs_comparation_res$specificity_mat[ln, col] <- res[2]
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
  row_names_numb = str_sub(row_names, 4, 7)
  col_names_numb = str_sub(col_names, 4, 7)
  
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
      i = which(str_sub(row_names,4,7) %in% bacs_matrix[bacs_matrix[, 2]==bac_i, 1])
      # On recupere ici les indices en ligne de toutes les billes faisant partie du bac j
      j = which(str_sub(col_names,4,7) %in% bacs_matrix[bacs_matrix[, 2]==bac_j, 1])
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

ncol = 562
ncells = 250
nrep = 40
data_dim = ncol*(ncol-1)/2
set.seed(99999)
nb_cluster = 250
resolution = "3Mb"
load("rdata/all_rda_data/umap.out_.rda")
load("rdata/all_rda_data/cellUperDiagData.rda")
load("rdata/all_rda_data/scHic_promoters_ids.rda")
load("rdata/all_rda_data/all_net_result_complex_3Mb_.rda")

all_net_result_complex_deg = all_net_result_complex_
for (bl in 2:16) {
  for (i in seq_len(length(all_net_result_complex_[[bl]]$crhs))) {
    mat = all_net_result_complex_[[bl]]$crhs[[i]]$mat_incidence
    if(sum(mat)!=-1){
      all_net_result_complex_deg[[bl]]$crhs[[i]][[2]]= degenerationMatrix(mat)
    }
  }
}


boot_fn = function(data, rows){
  print("Bootstrap beguin")
  cell_data = data[rows, ]
  umap_data = umap.out$layout[rows, ]
  # umap.out = umap(cell_data)
  
  cells_clusters <- kmeans(umap_data, nb_cluster, trace=FALSE)
  
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
  
  print("complete computing")
  
  return(c(median(highest_spec_), median(highest_sen_), mean(highest_spec_), mean(highest_sen_)))
}

r_boot = 100
ncpus = 4

tictoc::tic("Boot")
boot_res <- boot(cellUperDiagData, boot_fn, R=r_boot, parallel = "multicore", ncpus = ncpus)
tictoc::toc()

save(boot_res, file = "rdata/all_rda_data/boot_res.rda")

######################### Une estimation

cells_clusters <- kmeans(umap.out$layout, nb_cluster, trace=FALSE)

cluster_matrix_result = constr_mat_contacts(cells_clusters$cluster, cellUperDiagData)

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
summary(highest_spec_)

highest_sen = highest_sen[highest_sen[, 1] != "", ]
highest_sen_ = as.numeric(highest_sen[, 3])
summary(highest_sen_)

load("rdata/all_rda_data/boot_res.rda")

boot_res$t0

2*boot_res$t0 - colMeans(boot_res$t)
summary(boot_res$t[, 3])

boot.ci(boot_res, type=c("norm", "basic", "perc"), index = 3)
boot.ci(boot_res, type=c("norm", "basic", "perc"), index = 4)
plot(boot_res, index=4)

quantile(boot_res$t[, 4], probs = c(0.025, 0.975))

n = length(highest_spec_)
p = median(highest_spec_)
spec_IC_inf = p - 1.96*sqrt(p*(1-p)/n)
spec_IC_inf
spec_IC_sup = p + 1.96*sqrt(p*(1-p)/n)
spec_IC_sup

colMeans(boot_res$t)-boot_res$t0
apply(boot_res$t,2,sd)


# ###### Calcul des distances
# # Calcul de la distance avec ldes données de UMAP
# clusters_to_inspect = row.names(cellUperDiagData)[c(1:40, seq(40)+250, seq(40)+250*2, seq(40)+250*3, seq(40)+250*4, seq(40)+250*5, seq(40)+250*6, seq(40)+250*7, seq(40)+250*8, seq(40)+250*9, seq(40)+250*10)]
# cell_dist1 = as.matrix(dist(umap.out$layout[clusters_to_inspect, ], method = "euclidean"))
# cell_dist2 = as.matrix(dist(cellUperDiagData[clusters_to_inspect, ], method = "euclidean"))
# 
# pas = 40*39/2
# variances_intra_cells = c()
# variances_intra_cells_ = c()
# dt_variances_inter_cells = matrix(NA, nrow = pas*length(clusters_to_inspect)/40, ncol = 2)
# dt_variances_inter_cells_ = matrix(NA, nrow = pas*length(clusters_to_inspect)/40, ncol = 2)
# 
# for (clus in seq_len(length(clusters_to_inspect)/40)) {
#   cell_names = clusters_to_inspect[seq(40)*clus]
#   
#   mat = cell_dist1[
#     cell_names,
#     cell_names
#   ]
#   variances_intra_cells = c(
#     variances_intra_cells,
#     var(mat[upper.tri(mat)])
#   )
#   l = which(is.na(dt_variances_inter_cells[, 1]))[1]
#   
#   dt_variances_inter_cells[l:(l+pas-1), 1] = clus
#   dt_variances_inter_cells[l:(l+pas-1), 2] = mat[upper.tri(mat)]
#   
#   mat_ = cell_dist2[
#     cell_names,
#     cell_names
#   ]
#   variances_intra_cells_ = c(
#     variances_intra_cells_,
#     var(mat_[upper.tri(mat_)])
#   )
#   
#   dt_variances_inter_cells_[l:(l+pas-1), 1] = clus
#   dt_variances_inter_cells_[l:(l+pas-1), 2] = mat_[upper.tri(mat_)]
# }
# 
# grouped_data <- aggregate(dt_variances_inter_cells_[, 2], by = list(dt_variances_inter_cells_[, 1]), FUN = mean)
# variances_inter_cells_ <- var(grouped_data$x)
# 
# grouped_data <- aggregate(dt_variances_inter_cells[, 2], by = list(dt_variances_inter_cells[, 1]), FUN = mean)
# variances_inter_cells <- var(grouped_data$x)
# 
# variances_intra_cells_
# variances_inter_cells_
# 
# variances_intra_cells
# variances_inter_cells

load("rdata/all_rda_data/res_matrix.rda")
res = res_matrix



