<<<<<<< HEAD
source("./scripts/enhPromMousse6Mb.R")
source("./scripts/load_save_data.R")
source("./scripts/create_clusters.R")
source("./scripts/chevauche_CRHs_fn.R")
source("./scripts/crhs_candidates.R")
source("./scripts/create_3D_graphics_fn.R")

# source("./scripts/compute_distances_fn.R")
#source("./scripts/create_2D_graphics_fn.R")
# source("./scripts/chevauche_CRHs_fn.R")
# source("./scripts/find_min_by_rank.R")

# On peut voir que dans tous blocks, à des proportions différentes, on y trouve au moins un promoter
sort(unique(structure_1[promoters_ids, ]$X4))
table(structure_1[promoters_ids, ]$X4)

######################### Creation des CRHs #########################

### Structure 1 block 1
#structure_1_net <- create_clust_graph(structure_1, 1, 1:16, 3)
structure_1_net_bip <- create_bip_clust_graph(structure_1, promoters_ids, 1, 1:16, 3)
structure_1_ <- structure_1%>%
  left_join(
    structure_1_net$clust_block$structure,
    by = "ID"
  )%>%
  filter(!is.na(nb_edges))


### Structure 2 block 1
#structure_2_net <- create_clust_graph(structure_2, 2, 1:16, 3)
structure_2_net_bip <- create_bip_clust_graph(structure_2, promoters_ids, 2, 1:16, 3)
structure_2_ <- structure_2%>%
  left_join(
    structure_2_net$clust_block$structure,
    by = "ID"
  )%>%
  filter(!is.na(nb_edges))

################### Methodes d'appariement des CRHs ###################

### 1 - Similarité entre les CRHs en terme de nombre d'élements
overlap_crh <- chr_overlap(structure_1_net, structure_2_net)

### 2 - Similarité entre les CRHs en terme du nombre d'arrêtes des éléments
overlap_edge <- edge_identity_overlap(structure_1_net, structure_2_net)

### 3 - On récupère ici pour chaque block les CRHs qui été selectionnés par les deux méthodes 
crhs_select <- select_crh_candidats(overlap_crh, overlap_edge)

# Afficher les matchs selectionnés
for (i in 1:16) {
  print(crhs_select[[i]]$crh_edge_match)
}

### 4 Représentation graphique des CRHs matchés
### On choisi ici de representer les 4 premier dans les deux premiers replicats
grp <- plot_3D_plotly(structure_1_, structure_2_, 1, crhs_select)
grp







||||||| a6571cd
:diffg REsource("./scripts/load_save_data.R")
source("./scripts/create_clusters.R")
source("./scripts/chevauche_CRHs_fn.R")
source("./scripts/crhs_candidates.R")

source("./scripts/compute_distances_fn.R")
#source("./scripts/create_2D_graphics_fn.R")
source("./scripts/create_3D_graphics_fn.R")
# source("./scripts/chevauche_CRHs_fn.R")
source("./scripts/chevauche_CRHs_fn_.R")
source("./scripts/find_min_by_rank.R")

######################### Creation des CRHs #########################

### Structure 1 block 1
structure_1_net <- create_clust_graph(structure_1, 1, 1:16, 3)
structure_1_ <- structure_1%>%
  left_join(
    structure_1_net$clust_block$structure,
    by = "ID"
  )%>%
  filter(!is.na(nb_edges))


### Structure 2 block 1
structure_2_net <- create_clust_graph(structure_2, 2, 1:16, 3)
structure_2_ <- structure_2%>%
  left_join(
    structure_2_net$clust_block$structure,
    by = "ID"
  )%>%
  filter(!is.na(nb_edges))

################### Methodes d'appariement des CRHs ###################

<<<<<<< HEAD
### 1 - Similarité entre les CRHs en terme de nombre d'élements
overlap_crh <- chr_overlap(structure_1_net, structure_2_net)
=======
### 1- Méthode des distances
# Calcul des distances entre les CRHs
dist_r1_r2_bl1 <- m_dist_func(distance_diff1, distance_diff2)

# Distance minimum en ligne comme en colonne
dist_r1_r2_bl1_min <- cbind(
  sort(apply(dist_r1_r2_bl1, MARGIN = 1, FUN = min)),
  sort(apply(dist_r1_r2_bl1, MARGIN = 2, FUN = min))
)

dist_r1_r2_bl1_min

### 2 - Similarité entre les CRHs en terme de nombre d'élements
overlap <- chr_overlap(structure_1_1_net_comp, structure_2_1_net_comp)

overlap_perc <- 1 - overlap$chevauche_1_2_perc

overlap_perc[is.na(overlap_perc)] <- 99

overlap_min <- cbind(
  sort(apply(overlap_perc, MARGIN = 1, FUN = min, na.rm = TRUE)[apply(overlap_perc, MARGIN = 1, FUN = min, na.rm = TRUE) != 99]),
  sort(apply(overlap_perc, MARGIN = 2, FUN = min, na.rm = TRUE)[apply(overlap_perc, MARGIN = 2, FUN = min, na.rm = TRUE) != 99])
)

overlap_min
### 3 - Similarité entre les CRHs en terme du nombre d'arrêtes des éléments
overlap_edg <- edges_overlap(structure_1_1, structure_2_1, 1, 2, "nb_edges", "crh_id")

overlap_edg_max <- cbind(
  sort(apply(overlap_edg$chevauche, MARGIN = 1, FUN = max, na.rm = TRUE),decreasing = T),
  sort(apply(overlap_edg$chevauche, MARGIN = 2, FUN = max, na.rm = TRUE),decreasing = T)
)

overlap_edg <- edge_identity_overlap(structure_1_1_net_comp, structure_2_1_net_comp,structure_1_1_net$dist_bin, structure_2_1_net$dist_bin)

c1max = apply(overlap_edg$chevauche_1_2_perc,1,max,na.rm=T)
c1max[c1max==-Inf]=NA
c2max = apply(overlap_edg$chevauche_1_2_perc,2,max,na.rm=T)
c2max[c2max==-Inf]=NA

overlap_edg_max <- cbind(sort(c1max,decreasing = T),sort(c2max,decreasing = T)[1:sum(!is.na(c1max))])


overlap_edg_max
sum(overlap_edg_max[,1]==overlap_edg_max[,2])
### Comparaison des résulats des méthodes 
dist_r1_r2_bl1_min
overlap_min
overlap_edg_max


>>>>>>> f17f9cda862cba86b7b0c7ba5053125c8850617a

### 2 - Similarité entre les CRHs en terme du nombre d'arrêtes des éléments
overlap_edge <- edge_identity_overlap(structure_1_net, structure_2_net)


### 3 - On récupère ici pour chaque block les CRHs qui été selectionnés par les deux méthodes 
crhs_select <- select_crh_candidats(overlap_crh, overlap_edge)

# Afficher les matchs selectionnés
for (i in 1:16) {
  print(crhs_select[[i]]$crh_edge_match)
}

### 4 Représentation graphique des CRHs matchés
### On choisi ici de representer les 4 premier dans les deux premiers replicats
grp <- plot_3D_plotly(structure_1_, structure_2_, 1, crhs_select)
grp







=======
source("./scripts/enhProMousse.R")
source("./scripts/load_save_data.R")
source("./scripts/create_clusters.R")
source("./scripts/chevauche_CRHs_fn.R")
source("./scripts/crhs_candidates.R")
source("./scripts/create_3D_graphics_fn.R")

# source("./scripts/compute_distances_fn.R")
#source("./scripts/create_2D_graphics_fn.R")
# source("./scripts/chevauche_CRHs_fn.R")
# source("./scripts/find_min_by_rank.R")

# On peut voir que dans tous blocks, à des proportions différentes, on y trouve au moins un promoter
sort(unique(structure_1[promoters_ids, ]$X4))
table(structure_1[promoters_ids, ]$X4)

######################### Creation des CRHs #########################

### Structure 1 block 1
structure_1_net <- create_clust_graph(structure_1, 1, 1:16, 3)
structure_1_ <- structure_1%>%
  left_join(
    structure_1_net$clust_block$structure,
    by = "ID"
  )%>%
  filter(!is.na(nb_edges))


### Structure 2 block 1
structure_2_net <- create_clust_graph(structure_2, 2, 1:16, 3)
structure_2_ <- structure_2%>%
  left_join(
    structure_2_net$clust_block$structure,
    by = "ID"
  )%>%
  filter(!is.na(nb_edges))

################### Methodes d'appariement des CRHs ###################

### 1 - Similarité entre les CRHs en terme de nombre d'élements
overlap_crh <- chr_overlap(structure_1_net, structure_2_net)

### 2 - Similarité entre les CRHs en terme du nombre d'arrêtes des éléments
overlap_edge <- edge_identity_overlap(structure_1_net, structure_2_net)

### 3 - On récupère ici pour chaque block les CRHs qui été selectionnés par les deux méthodes 
crhs_select <- select_crh_candidats(overlap_crh, overlap_edge)

# Afficher les matchs selectionnés
for (i in 1:16) {
  print(crhs_select[[i]]$crh_edge_match)
}

### 4 Représentation graphique des CRHs matchés
### On choisi ici de representer les 4 premier dans les deux premiers replicats
grp <- plot_3D_plotly(structure_1_, structure_2_, 1, crhs_select)
grp







>>>>>>> c68847982f270febe8e97d841a91df2dfde1aad9
