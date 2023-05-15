# Chargement des librairies
source("scripts/crhs_creation_with_cells/cells_clustering")
source("scripts/crhs_creation_with_cells/matrix_for_juicer_creation.R")


# Chargement des librairies
library(igraph)
library(dplyr)
library(stringr)
library(stats)
library(data.table)
library(tidyverse)

# Chargement des promoters
load("rdata/scHic_promoters_ids.rda")

# Chargement des données sur le clustering
load("rdata/cellUperDiagData.rda")
cellUperDiagData = t(cellUperDiagData)
row.names(cellUperDiagData) <- str_c("cellule_", 1:nrow(cellUperDiagData))

save(clu_chrs, file = "rdata/clu_chrs_result.rda")

