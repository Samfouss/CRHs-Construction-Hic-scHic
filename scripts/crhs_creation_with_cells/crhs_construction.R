
# Chargement des données sur les clusters
load("rdata/all_rda_data/cluster_matrix_result.rda")
source("./scripts/crhs_creation_with_cells/create_graph_from_cells_fn.R")

clus_num = c()
for (cell in names(cluster_matrix_result)) {
  n = as.numeric(unlist(strsplit(cell, "[_]"))[2])
  clus_num = c(clus_num, n)
}

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
  
  names(clu_chrs_result) <- str_c("cluster_", clus_num)
  clu_chrs_result
}

resolution = "3Mb"

clu_chrs_result <- get_clusters_crhs(cluster_matrix_result, resolution = resolution)

############## Distribution des crhs ##########
ncrhs <- matrix(
  0,
  nrow = length(clu_chrs_result),
  ncol = 1,
)
for (res in seq_len(length(clu_chrs_result))) {
  # sum(fireCaller_result[[res]]$FIRE_output$Mousse_cells_clus_12_indicator)
  ncrhs[res, 1] <- length(clu_chrs_result[[res]])
}
sum(ncrhs[, 1])
summary(ncrhs[, 1])

ncells = 250
crhs_inspection = matrix(
  0,
  nrow = sum(ncrhs[, 1]),
  ncol = 3
)

l = 1

for (cell in seq_len(ncells)) {
  for (j in seq_len(length(clu_chrs_result[[cell]]))) {
    crhs_inspection[l, 1] =  dim(clu_chrs_result[[cell]][[j]]$mat_incidence)[1]
    crhs_inspection[l, 2] =  dim(clu_chrs_result[[cell]][[j]]$mat_incidence)[2]
    crhs_inspection[l, 3] =  paste0("cell ", sprintf("%03d", clus_num[cell]), " - CRH ", j)
    l = l +  1
  }
}

plot(
  table(as.numeric(crhs_inspection[, 1])),
  xlab = "Nombre de CRHs dans les structures",
  ylab = "Nombre de répétition"
)
plot(
  table(as.numeric(crhs_inspection[, 2])),
  xlab = "Nombre de CRHs dans les structures",
  ylab = "Nombre de promoters"
)

## ================================================ Retenir les CRHs complexes ===========
clu_chrs_result_complex = clu_chrs_result
for (cell in seq_len(ncells)) {
  for (j in seq_len(length(clu_chrs_result[[cell]]))) {
    mat = clu_chrs_result_complex[[cell]][[j]]$mat_incidence
    
    # Identification des CRHs moins complex
    if((nrow(mat) == 1 & ncol(mat) == 1) | (nrow(mat) == 2 & ncol(mat) == 1) | (nrow(mat) == 1 & ncol(mat) == 2)){
      clu_chrs_result_complex[[cell]][[j]]$mat_incidence <- -1
    }
  }
}

crhs_inspection_ = matrix(
  0,
  nrow = sum(ncrhs[, 1]),
  ncol = 3
)

l = 1

for (cell in seq_len(ncells)) {
  for (j in seq_len(length(clu_chrs_result_complex[[cell]]))) {
    mat = clu_chrs_result_complex[[cell]][[j]]$mat_incidence
    if(sum(mat)!=-1){
      crhs_inspection_[l, 1] =  dim(mat)[1]
      crhs_inspection_[l, 2] =  dim(mat)[2]
      crhs_inspection_[l, 3] =  paste0("cell ", sprintf("%03d", clus_num[cell]), " - CRH ", j)
      l = l +  1
    }
  }
}

crhs_inspection_ = crhs_inspection_[crhs_inspection_[, 1]!="0", ]

plot(
  table(as.numeric(crhs_inspection_[, 1])),
  xlab = "Nombre de CRHs dans les structures",
  ylab = "Nombre de répétition"
)
plot(
  table(as.numeric(crhs_inspection_[, 2])),
  xlab = "Nombre de CRHs dans les structures",
  ylab = "Nombre de promoters"
)
## ================================================ Retenir les CRHs complexes ===========

########### Sauvegarde des données ###########
if(resolution=="2Mb"){
  save(clu_chrs_result, file = "rdata/all_rda_data/clu_chrs_result_2Mb.rda")
}else if(resolution=="3Mb"){
  save(clu_chrs_result, file = "rdata/all_rda_data/clu_chrs_result_3Mb.rda")
}else{
  save(clu_chrs_result, file = "rdata/all_rda_data/clu_chrs_result_6Mb.rda")
}
