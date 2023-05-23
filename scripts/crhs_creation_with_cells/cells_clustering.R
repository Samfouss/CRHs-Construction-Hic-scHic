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
row.names(cellUperDiagData) <- str_c("cellule_", 1:nrow(cellUperDiagData))

load("rdata/cellUperDiagData_with_rep.rda")
cellUperDiagData_with_rep = t(cellUperDiagData_with_rep)
row.names(cellUperDiagData_with_rep) <- str_c("cellule_", 1:nrow(cellUperDiagData_with_rep))

############################ Kmeans avec 5 class ######################

# Afin de reproduire les mêmes résulats en choisissant les memes points de départ (puisqu'ils sont aléatoire)
set.seed(99999)
nb_class = 4
cells_clusters_4 <- kmeans(
  cellUperDiagData_with_rep,
  nb_class,
  algorithm = c("Hartigan-Wong"),
  iter.max = 50,
  trace=FALSE
)
cells_clusters_4$withinss
cells_clusters_4$betweenss

constr_mat_contacts(cells_clusters_4$cluster, cellUperDiagData_with_rep)

############################ Kmeans avec 4 class ######################

nb_class = 5
cells_clusters_5 <- kmeans(
  cellUperDiagData_with_rep,
  nb_class,
  algorithm = c("Hartigan-Wong"),
  iter.max = 50,
  trace=FALSE
)
cells_clusters_5$withinss
cells_clusters_5$betweenss

constr_mat_contacts(cells_clusters_5$cluster, cellUperDiagData_with_rep)

############################ Kmeans avec 3 class ######################

nb_class = 3
cells_clusters <- kmeans(
  cellUperDiagData,
  nb_class,
  algorithm = c("Hartigan-Wong"),
  iter.max = 50,
  trace=FALSE
)
cells_clusters$withinss
cells_clusters$betweenss


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