
# Chargement des données sur les clusters
load(paste0("rdata/all_rda_data/merge_loops_clus_", nb_clusters, ".rda"))
load("rdata/all_rda_data/cluster_matrix_result.rda")

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
  
  names(clu_chrs_result) <- str_c("cluster", seq_len(length(clusters_matrix)))
  clu_chrs_result
}

clu_chrs_result <- get_clusters_crhs(cluster_matrix_result, resolution = "2Mb")

############## Distribution des crhs ##########
ncrhs <- matrix(
  0,
  nrow = length(clu_chrs_result),
  ncol = 1,
)
for (res in seq_len(length(clu_chrs_result))) {
  # sum(fireCaller_result[[res]]$FIRE_output$Mousse_cells_clus_12_indicator)
  n = 0
  for (i in seq_len(length(clu_chrs_result[[res]]))) {
    n = n + length(clu_chrs_result[[res]][[i]]$crhs)
  }
  ncrhs[res, 1] <- n
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
    for (i in seq_len(length(clu_chrs_result[[cell]][[j]]$crhs))) {
      crhs_inspection[l, 1] =  dim(clu_chrs_result[[cell]][[j]]$crhs[[i]]$mat_incidence)[1]
      crhs_inspection[l, 2] =  dim(clu_chrs_result[[cell]][[j]]$crhs[[i]]$mat_incidence)[2]
      crhs_inspection[l, 3] =  paste0("cell ", sprintf("%03d", cell), "- CRH ", i)
      l = l +  1
    }
  }
}

########### Sauvegarde des données ###########
save(clu_chrs_result, file = "rdata/all_rda_data/clu_chrs_result.rda")

