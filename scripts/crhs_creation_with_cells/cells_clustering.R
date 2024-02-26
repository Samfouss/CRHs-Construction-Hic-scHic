library(ggpubr)
library(factoextra)



source("scripts/crhs_creation_with_cells/juicerInputFileCreation_fn.R")

# Chargement des promoters
load("rdata/all_rda_data/scHic_promoters_ids.rda")

# dataToUse = rep_data; ideal_rep_data; HicImputeDtata
# clustering_meth = kmeans; KMedoide; ACP
# nb_cluster

cells_lustering <- function(dataToUse = "ideal_rep_data", clustering_meth = "KMedoide", kmeans_algo = "Hartigan-Wong" , nb_cluster = 6){
  
  # Chargement des données sur le clustering
  if(dataToUse=="rep_data"){
    
    load("rdata/all_rda_data/cellUperDiagData_with_repT.rda")
    ncell = 250
    cellUperDiagData = t(cellUperDiagData_with_repT)
    row.names(cellUperDiagData) <- paste0("cellule_", seq_len(nrow(cellUperDiagData)), "_rep", seq_len(ncell))
    
  }else if(dataToUse=="ideal_rep_data"){
    
    load("rdata/all_rda_data/cellUperDiagData_with_rep.rda")
    cellUperDiagData = t(cellUperDiagData_with_rep)
    row.names(cellUperDiagData) <- str_c("cellule_", 1:nrow(cellUperDiagData))
    
  }else if(dataToUse=="HicImputeDtata"){
    
    load("rdata/all_rda_data/MCMCImpute_result.rda")
    cellUperDiagData = t(MCMCImpute_result$Impute_SZ)
    row.names(cellUperDiagData) <- str_c("cellule_", 1:nrow(cellUperDiagData))
    
  }
  
  if(clustering_meth == "kmeans"){
    
    set.seed(99999)
    cells_clusters <- kmeans(
      cellUperDiagData,
      nb_cluster,
      algorithm = kmeans_algo,
      iter.max = 100,
      trace=FALSE
    )
    print(cells_clusters$size)
    constr_mat_contacts(cells_clusters$cluster, cellUperDiagData)
    
  }else if(clustering_meth == "KMedoide"){
    
    cells_clusters = pamk(cellUperDiagData, nb_cluster)
    table(cells_clusters$clustering)
    constr_mat_contacts(cells_clusters$pamobject$clustering, cellUperDiagData)
    
  }else if(clustering_meth == "ACP"){
    
    nb_class_max = nb_cluster
    nb_class_min = 5
    res.pca <- PCA(cellUperDiagData, graph = FALSE)
    # eig.val <- get_eigenvalue(res.pca)
    # fviz_mca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, ggtheme = theme_minimal())
    ncp_choose = which(res.pca$eig[, 3, drop = FALSE]>=80)[3]
    res.pca <- PCA(cellUperDiagData, ncp = ncp_choose, graph = FALSE)
    cells_clusters <- HCPC(res.pca, kk = Inf, min = nb_class_min, max = nb_class_max, graph = FALSE, consol = TRUE)
    # fviz_cluster(cells_clusters, repel = TRUE, geom = "point", main = "Classification")
    
    # Nombre de cellules par cluster
    table(cells_clusters$data.clust$clust)
    clusters = unique(cells_clusters$data.clust$clust)
    
    # Résumé
    summary(as.vector(table(cells_clusters$data.clust[, ncol(cells_clusters$data.clust), drop = FALSE])))
    
  }
  
  save(cells_clusters, file = "rdata/all_rda_data/cluster_result.rda")
  
  cells_clusters
}




################## clustering hierarchique #####################

# hclust_clusters = hclust(cellUperDiagData)

################## dbscan clustering #####################

# dbscan_clusters = dbscan(cellUperDiagData, eps = 0.2, MinPts = 10)


#################################### Distribution ####################
# nb_class = 6
# for (clus in seq_len(nb_class)) {
#   print(paste0("cluster ", clus))
#   #dtl = read.table(paste0("../BackUp14022020/ClusterOutput/Cluster_6/input/cluster_", clus, "_juicer_input.txt"))
#   dtl = read.table(paste0("rdata/juicerInputFiles/cluster_", clus, "_juicer_input.txt"))
#   print(dim(dtl))
#   print(summary(dtl[, 9]))
# }




