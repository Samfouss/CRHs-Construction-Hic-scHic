# Exécution des scripts nécéssaires
source("./scripts/load_save_data.R")
source("./scripts/create_clusters.R")
source("./scripts/chevauche_CRHs_fn.R")
source("./scripts/fusionne_replicats_fn.R")

# On peut voir que dans tous blocks, à des proportions différentes, on y trouve au moins un promoter
sort(unique(structure_1[promoters_ids, ]$X4))
promoters_ids[which(promoters_ids<341)]

# Un petit resumé du nombre de promoters par bloc
table(structure_1[promoters_ids, ]$X4)
table(structure_2[promoters_ids, ]$X4)
table(structure_3[promoters_ids, ]$X4)

######################### Creation des CRHs #########################

# Test de fusion de trois replicats

net1 <- create_bip_clust_graph2(structure_1, promoters_ids, 1, 1:16, 3)
net2 <- create_bip_clust_graph2(structure_2, promoters_ids, 2, 1:16, 3)
net3 <- create_bip_clust_graph2(structure_3, promoters_ids, 3, 1:16, 3)

edgover <- edge_identity_overlap2(net1, net2)

fusion_essaie <- fusion_comp2(net1, net2, edgover)

# On a ici le résumé des fusions
# Chaque ligne correspond aux CRHs fusionnés et non fusionnés
# - La première valeure est le numero du replica
# - La deuxième valeure est le numero du CRHs dans le replica
fusion_essaie$block1$resume_fusion

edgover <- edge_identity_overlap2(net3, fusion_essaie)

fusion_essaie <- fusion_comp2(net3, fusion_essaie, edgover)






# On récupère ici le nombre de réplicat
nb_replicas = length(unique(all_single_structure$replica))

# Bien avant de rentrer dans la boucle, on initialise la matrice de fusion avec le cluster de la première structure
all_net_result <- create_bip_clust_graph2(
  as_tibble(
    all_single_structure%>%
      filter(replica == sprintf("%03d", 1))%>%
      select(-c(ends_with("_c"), "replica"))%>%
      mutate(
        ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n()))
      )
  ),
  promoters_ids, 
  3, 
  1:16, 
  3
)

nb_replicas = 50
for (r in 2:nb_replicas) {
  structure <- as_tibble(
    all_single_structure%>%
      filter(replica == sprintf("%03d", r))%>%
      select(-c(ends_with("_c"), "replica"))%>%
      mutate(
        ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n()))
      )
  )
  
  print(r)
  print("Cluster")
  net <- create_bip_clust_graph2(structure, promoters_ids, r, 1:16, 3)
  
  print("Overlap")
  overlap_edge <- edge_identity_overlap2(
    net, 
    all_net_result
  )
  
  print("Fusion")
  all_net_result <- fusion_comp2(net, all_net_result, overlap_edge)
  
}




