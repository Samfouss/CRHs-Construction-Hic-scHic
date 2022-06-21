source("./scripts/enhPromMousse6Mb.R")
source("./scripts/load_save_data.R")
source("./scripts/create_clusters.R")
source("./scripts/chevauche_CRHs_fn.R")
source("./scripts/crhs_candidates.R")

source("./scripts/compute_distances_fn.R")
#source("./scripts/create_2D_graphics_fn.R")
source("./scripts/create_3D_graphics_fn.R")
# source("./scripts/chevauche_CRHs_fn.R")
source("./scripts/chevauche_CRHs_fn_.R")
source("./scripts/find_min_by_rank.R")

# source("./scripts/compute_distances_fn.R")
#source("./scripts/create_2D_graphics_fn.R")
# source("./scripts/chevauche_CRHs_fn.R")
# source("./scripts/find_min_by_rank.R")

# On peut voir que dans tous blocks, à des proportions différentes, on y trouve au moins un promoter
sort(unique(structure_1[promoters_ids, ]$X4))
table(structure_1[promoters_ids, ]$X4)

######################### Creation des CRHs #########################

### Structure 1 block 1
# Création des CRHs avec la structure 1
structure_1_net_bip <- create_bip_clust_graph(structure_1, promoters_ids, 1, 1:16, 3)

### Structure 2 block 1
# Création des CRHs avec la structure 2
structure_2_net_bip <- create_bip_clust_graph(structure_2, promoters_ids, 2, 1:16, 3)

################### Methodes d'appariement des CRHs ###################

### Similarité entre les CRHs en terme du nombre d'arrêtes des éléments
overlap_edge <- edge_identity_overlap(structure_1_net_bip, structure_2_net_bip)

# Afficher les matchs selectionnés
for (i in 1:16) {
  print(paste0("#### Les CRHs qui matchent dans les deux sens (en ligne ou en colonne) dans le block ", i, " ####"))
  print(overlap_edge[[i]]$chev_edge_comp)
}

### 4 Représentation graphique des CRHs matchés
### On choisi ici de representer les 4 premier dans les deux premiers replicats
grp <- plot_3D_plotly(structure_1_, structure_2_, 1, crhs_select)
grp

