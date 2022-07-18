# Exécution des scripts nécéssaires
source("./scripts/load_save_data.R")
source("./scripts/create_clusters.R")
source("./scripts/chevauche_CRHs_fn.R")

# On peut voir que dans tous blocks, à des proportions différentes, on y trouve au moins un promoter
sort(unique(structure_1[promoters_ids, ]$X4))
promoters_ids[which(promoters_ids<341)]

# Un petit resumé du nombre de promoters par bloc
table(structure_1[promoters_ids, ]$X4)
table(structure_2[promoters_ids, ]$X4)
table(structure_3[promoters_ids, ]$X4)
table(structure_4[promoters_ids, ]$X4)

######################### Creation des CRHs #########################

# Test de fusion de trois replicats

net1 <- create_bip_clust_graph(structure_1, promoters_ids, 1, 1:16, 3)
net2 <- create_bip_clust_graph(structure_2, promoters_ids, 2, 1:16, 3)
net3 <- create_bip_clust_graph(structure_3, promoters_ids, 3, 1:16, 3)


overlap_edge <- edge_identity_overlap(
  net1, 
  net2
)

result <- fusion_comp(net1, net2, overlap_edge)

