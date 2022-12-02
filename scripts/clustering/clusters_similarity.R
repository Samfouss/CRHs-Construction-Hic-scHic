# Chargement des librairies
library("devtools")
## Warning: le package 'devtools' a été compilé avec la version R 4.1.3
## Le chargement a nécessité le package : usethis
## Warning: le package 'usethis' a été compilé avec la version R 4.1.3
# Install "HiCImpute" package from github.
# install_github("https://github.com/sl-lin/HiCImpute")
library("HiCImpute")

library("tidyverse")
library("FactoMineR")
library("factoextra")
library("ggpubr")

nb_class_max = 50
nb_class_min = 5


# Chargement des données sur les cellules et construction du clustering de cellules
load("rdata/cells_matrix.rda")
cells_data <- cells_matrix[[1]]
for (b in 2:16) {
  cells_data <- cbind(
    cells_data,
    cells_matrix[[b]] 
  )
}

cells_data[cells_data == 1] <- "in"
cells_data[cells_data == 0] <- "Not in"

cells_data <- data.frame(cells_data)

res.mca <- MCA(cells_data, graph = FALSE)
eig.val <- get_eigenvalue(res.mca)
fviz_mca_ind(res.mca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, ggtheme = theme_minimal())
which(res.mca$eig[, 3, drop = FALSE]>=70)[1:5]

res.mca <- MCA(cells_data, ncp = 249, graph = FALSE)
hcpc_cluster <- HCPC(res.mca, kk=Inf, min = nb_class_min, max = nb_class_max, graph = FALSE)
fviz_cluster(hcpc_cluster, repel = TRUE, geom = "point", main = "Classification des cellules après avoir retiré \n les CRHs moins complexes")

table(hcpc_cluster$data.clust$clust)


# Disposition des données dans une liste pour le traitement avec HicImpute
ncol = 562
ncells = 250
data_dim = ncol*(ncol-1)/2

hicImpute_data <- matrix(
  0, 
  nrow = data_dim, 
  ncol = ncells, 
  byrow = TRUE
)

for (i in seq_len(ncells)) {
  mat <- as.matrix(
    read.table(
      paste0("rdata/single_cell_hic_data/hic_mat_", sprintf("%03d", i), ".txt"), 
      quote="\"", 
      comment.char="", 
      stringsAsFactors = FALSE
    ) 
  )
  
  hicImpute_data[, i] <- mat[upper.tri(mat)]
}
rm("mat")
colnames(hicImpute_data) <- str_c("cell", sprintf("%03d", seq(ncells)))


# Nous prenons ici les paramètres par défaut de la fonction MCMCImpute. Cepandant comme nous ne disposons pas de données en vrac, nous mettons bulk à NULL. Nous pouvons cependant aussi générer les données en vrac à partir du code du simulation de départ en prennant en compte toutes les cellules comme le faisait les auteurs.
set.seed(1234)

max_cluters_hicImput = 42
cluster=scHiC_Kmeans(
  hicImpute_data, 
  centers=max_cluters_hicImput, 
  nstart=1, 
  iter.max=1000, 
  seed=1
)

table(cluster$cluster)

cluster$cluster[cluster$cluster==1]

# Dans cluster_data, on résume les données sur les clustering
cluster_data <- as_tibble(
  data.frame(
    hicImput_clus = cluster$cluster,
    hicImput_cells = names(cluster$cluster),
    hcpc_clus = hcpc_cluster$data.clust$clust,
    hcpc_cells = row.names(hcpc_cluster$data.clust)
  )
)


hicImput_cluster = unique(cluster_data$hicImput_clus)
hcpc_cluster = unique(cluster_data$hcpc_clus)

similarity_results = matrix(
  NA,
  nrow = length(hicImput_cluster), 
  ncol = 5
)

# Dans cette boucle, on parcour les clusters crées par HCPC et HicImput puis calcul la similarité entre les clusters. L'objectif est de voir dans les deux cas, quels sont les clusters qui sont similaires
for (i in seq(length(hcpc_cluster))) {
  
  # k parcours les différents clusters
  k = hcpc_cluster[i]
  # cellsi recupère les elements du cluster
  cellsi = cluster_data[cluster_data$hcpc_clus==k, ]$hcpc_cells
  # perc : est initialisé à 0 pour à chaque fois recupérer le % de similarité entre les clusters
  perc = 0
  for (j in hicImput_cluster) {
    
    cellsj = cluster_data[cluster_data$hicImput_clus==j, ]$hicImput_cells
    perc = length(intersect(cellsi, cellsj))/max(length(cellsi), length(cellsj))
    
    if(perc>0){
      cluster_to_keep = c(j, length(cellsj), perc)
    }
  }
  
  similarity_results[i, 1] <- i
  similarity_results[i, 2] <- length(cellsi)
  similarity_results[i, 3] <- cluster_to_keep[1]
  similarity_results[i, 4] <- cluster_to_keep[2]
  similarity_results[i, 5] <- cluster_to_keep[3]
  
  hicImput_cluster = hicImput_cluster[hicImput_cluster != cluster_to_keep[1]]
}

similarity_results