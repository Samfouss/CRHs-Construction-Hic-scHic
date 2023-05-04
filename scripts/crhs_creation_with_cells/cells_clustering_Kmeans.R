library(stats)
library(stringr)

# Chargement des données sur le clustering
load("rdata/cellUperDiagData.rda")
cellUperDiagData = t(cellUperDiagData)
row.names(cellUperDiagData) <- str_c("cellule_", 1:nrow(cellUperDiagData))

# Afin de reproduire les mêmes résulats en choisissant les memes points de départ (puisqu'ils sont aléatoire)
set.seed(99999)
nb_class = 8
cells_clusters <- kmeans(
  cellUperDiagData,
  nb_class,
  algorithm = c("Hartigan-Wong"),
  iter.max = 50,
  trace=FALSE
)
cells_clusters$withinss
cells_clusters$betweenss