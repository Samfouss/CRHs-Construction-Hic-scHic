# Chargement des librairies
library(igraph)
library(dplyr)
library(stringr)
library(stats)
library(data.table)
library(tidyverse)
library(IRanges)
library(readr)

# Etape 1 : Transformation des données des matrices scHic
# source("scripts/crhs_creation_with_cells/getMatrixUperDataAndSave.R")

# Etape 2 :  Création de clusters
source("scripts/crhs_creation_with_cells/cells_clustering.R")
res = cells_lustering(dataToUse = "ideal_rep_data", clustering_meth = "KMedoide", nb_cluster = nb_clusters)

# Etape 3.1 : Transformation des sorties de Juicer en matrice d'incidence et sauvegarde des résultats
source("scripts/crhs_creation_with_cells/fromJuicerOutputFileToMatrix.R")

# Etaoe 3.2 : Transformation des sorties du clustering en matric d'incidence directement et sauvegarde des résualts
source("scripts/crhs_creation_with_cells/fromClusterToMatrix.R")

# Etape 4 : Chargement de la fonction pour la création des graphs
source("scripts/crhs_creation_with_cells/create_graph_from_cells_fn.R")

# Construction des réseaux et sauvegarde des résultats
source("scripts/crhs_creation_with_cells/crhs_construction.R")


