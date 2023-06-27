library(stats)
library(stringr)
library("FactoMineR")
library("factoextra")
source("scripts/crhs_creation_with_cells/juicerInputFileCreation_fn.R")

# Chargement des promoters
load("rdata/scHic_promoters_ids.rda")

# Chargement des données sur le clustering
load("rdata/cellUperDiagData.rda")
cellUperDiagData = t(cellUperDiagData)

load("rdata/cellUperDiagData_with_rep.rda")
cellUperDiagData = t(cellUperDiagData_with_rep)

load("rdata/MCMCImpute_result.rda")
cellUperDiagData = t(MCMCImpute_result$Impute_SZ)


row.names(cellUperDiagData) <- str_c("cellule_", 1:nrow(cellUperDiagData))
############################ Kmeans classes balancées ######################
start_cls = 5
end_cls = 15
set.seed(99999)
for (cls in seq(start_cls, end_cls)) {
  cells_clusters <- kmeans(
    cellUperDiagData,
    cls,
    algorithm = c("Hartigan-Wong"),
    iter.max = 50,
    trace=FALSE
  )
  
  print(paste0("Repartition des cellules dans les classes avec ", cls, " clusters"))
  print(cells_clusters$size)
}
# Afin de reproduire les mêmes résulats en choisissant les memes points de départ (puisqu'ils sont aléatoire)
nb_class = 10
set.seed(99999)
cells_clusters <- kmeans(
  cellUperDiagData,
  nb_class,
  algorithm = c("Hartigan-Wong"),
  iter.max = 50,
  trace=FALSE
)
cells_clusters$size
cells_clusters$withinss
cells_clusters$betweenss

points(cells_clusters$centers)

constr_mat_contacts(cells_clusters$cluster, cellUperDiagData)

################## K Medoids clustering #####################
library(fpc)

pamk_clusters = pamk(cellUperDiagData)
table(pamk_clusters$pamobject$clustering)

library(cluster)
pamk_clusters_8 = pam(cellUperDiagData, 8)
table(pamk_clusters_8$clustering)
constr_mat_contacts(pamk_clusters_8$clustering, cellUperDiagData)


pamk_clusters_5 = pamk(cellUperDiagData, 5)
table(pamk_clusters_5$pamobject$clustering)
constr_mat_contacts(pamk_clusters_5$pamobject$clustering, cellUperDiagData)

################## clustering hierarchique #####################
library(stats)
hclust_clusters = hclust(cellUperDiagData)

################## clustering hierarchique #####################

dbscan_clusters = dbscan(cellUperDiagData, eps = 0.2, MinPts = 10)


#################################### Distribution ####################

for (clus in seq_len(nb_class)) {
  print(paste0("cluster ", clus))
  dtl = read.table(paste0("rdata/juicerInputFiles/cluster_", clus, "_juicer_input.txt"))
  print(dim(dtl))
  print(summary(dtl[, 9]))
}

############################ Utilisation de l'ACP  ######################

res.pca <- PCA(cellUperDiagData, graph = FALSE)
# eig.val <- get_eigenvalue(res.pca)
# fviz_mca_ind(res.pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, ggtheme = theme_minimal())
ncp_choose = which(res.pca$eig[, 3, drop = FALSE]>=80)[3]
res.pca <- PCA(cellUperDiagData, ncp = ncp_choose, graph = FALSE)

nb_class_max = 8
nb_class_min = 3
hcpc_cluster <- HCPC(res.pca, kk = Inf, min = nb_class_min, max = nb_class_max, graph = FALSE, consol = TRUE)
fviz_cluster(hcpc_cluster, repel = TRUE, geom = "point", main = "Classification")

# Nombre de cellules par cluster
table(hcpc_cluster$data.clust$clust)
clusters = unique(hcpc_cluster$data.clust$clust)

# Résumé
summary(as.vector(table(hcpc_cluster$data.clust[, ncol(hcpc_cluster$data.clust), drop = FALSE])))











