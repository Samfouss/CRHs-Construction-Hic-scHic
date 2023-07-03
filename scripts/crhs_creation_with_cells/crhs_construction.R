
# Chargement des données sur les clusters
load(paste0("rdata/all_rda_data/merge_loops_clus_", nb_clusters, ".rda"))
load("rdata/all_rda_data/scHic_promoters_ids.rda")

clu_chrs_result = list()
for (clus in seq_len(length(merge_loops_clus))) {
  net = create_bip_clust_graph_from_cell(merge_loops_clus[[clus]], scHic_promoters_ids, clus)
  clu_chrs_result[[length(clu_chrs_result)+1]] <- net
}

names(clu_chrs_result) <- str_c("cluster", seq_len(length(merge_loops_clus)))

############## Distribution des crhs ##########
ncrhs <- matrix(
  0,
  nrow = length(clu_chrs_result),
  ncol = 1,
)
for (res in seq_len(length(clu_chrs_result))) {
  # sum(fireCaller_result[[res]]$FIRE_output$Mousse_cells_clus_12_indicator)
  ncrhs[res, 1] <- length(clu_chrs_result[[res]]$crhs)
}

summary(ncrhs[, 1])

########### Sauvegarde des données ###########
save(clu_chrs_result, file = paste0("rdata/all_rda_data/clu_chrs_result_", length(merge_loops_clus), ".rda"))

