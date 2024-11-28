set.seed(123)

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

library("tidyverse")
library("igraph")
library("boot")

# source("crhs_processing.R")
load("rdata/all_rda_data/all_net_result_complex_3Mb_.rda")
load("rdata/all_rda_data/cellUperDiagData.rda")
load("rdata/all_rda_data/scHic_promoters_ids.rda")

resolution = "3Mb"
num_groups = 250
group_size = 40
ncol = 562
ncells = 250
nrep = 40
data_dim = ncol*(ncol-1)/2
nb_cluster = 250
resolution = "3Mb"
cellUperDiagData_new <- matrix(NA, nrow = nrow(cellUperDiagData), ncol = ncol(cellUperDiagData))
r_boot = 1000
ncpus = 31

# Etape 3 : fonction permettant de faire la réduction de matrices
source("scripts/crhs_comparaison/degenerationMatrix_fn.R")

# Etape 4 : fonction permettant de faire la comparaison de CRHs
source("scripts/crhs_comparaison/compute_comparaison.R")

source("scripts/crhs_creation_with_cells/create_graph_from_cells_fn.R")

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
  print(paste0("Début : ", format(Sys.time(), "%H:%M:%S")))
  j = 1
  for (k in 1:num_groups) {
    # Définir les groupe de type de cellules : On a exactement 40 cellules par types de cellules
    row_range <- ((k - 1) * group_size + 1):(k * group_size)
    # Récupérer les cellules qui ont été tirées et le groupe auquel ils appartiennent
    common_numbers <- intersect(rows, row_range)
    # On somme par ligne les valeurs des matrices des cellules appartenant à la même cellule
    if(length(common_numbers)>0){
      cellUperDiagData_new[j, 1:ncol(cellUperDiagData_new)] <- colSums(cellUperDiagData[common_numbers, ])
      j = j + 1
    }
  }
  
  cellUperDiagData_new <- cellUperDiagData_new[complete.cases(cellUperDiagData_new), ]
  cellUperDiagData_new <- ifelse(cellUperDiagData_new>0, 1, 0)
  
  clusters_matrix <- list()
  for (cell in seq_len(num_groups)) {
    
    mat_row_name = paste0("BAC", sprintf("%03d", seq_len(ncol)))
    mat_temp <- matrix(
      0,
      nrow = length(mat_row_name),
      ncol = length(mat_row_name),
      dimnames = list(mat_row_name, mat_row_name)
    )
    
    mat_temp[upper.tri(mat_temp, diag = FALSE)] <- cellUperDiagData_new[cell, , drop = FALSE]
    mat_temp <- mat_temp + t(mat_temp)
    
    # mat_temp <- apply(mat_temp, 2, function(x) as.numeric(x>0))
    clusters_matrix[[length(clusters_matrix) + 1]] <- mat_temp
  }
  
  clu_chrs_result <- get_clusters_crhs(clusters_matrix, resolution = resolution)
  
  crhs_comparaison_res = compute_comparaison(all_net_result_complex_deg, clu_chrs_result)
  
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
  
  
  return(c(median(highest_spec_), median(highest_sen_), mean(highest_spec_), mean(highest_sen_)))
  print(paste0("Fin : ", format(Sys.time(), "%H:%M:%S")))
}


boot_res <- boot(cellUperDiagData, boot_fn, R=2)


boot_res <- boot(cellUperDiagData, boot_fn, R=r_boot, parallel = "multicore", ncpus = ncpus)

save(boot_res, file = "rdata/ideal_boot_res.rda")
print("End of bootstrap computation")



load("rdata/all_rda_data/ideal_boot_res.rda")
2*boot_res$t0 - colMeans(boot_res$t)
boot.ci(boot_res, type=c("norm", "basic", "perc"), index = 3)
boot.ci(boot_res, type=c("norm", "basic", "perc"), index = 4)



