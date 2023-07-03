# Exécution des scripts nécéssaires
source("./scripts/load_save_data.R")
source("./scripts/crhs_creation_with_bins/create_bip_graphs_fn.R")
source("./scripts/crhs_creation_with_bins/chevauche_CRHs_fn.R")
source("./scripts/crhs_creation_with_bins/fusionne_replicats_fn.R")

# On peut voir que dans tous blocks, à des proportions différentes, on y trouve au moins un promoter
sort(unique(structure_1[promoters_ids, ]$X4))
promoters_ids[which(promoters_ids<341)]

# Un petit resumé du nombre de promoters par bloc
table(structure_1[promoters_ids, ]$X4)
table(structure_2[promoters_ids, ]$X4)
table(structure_3[promoters_ids, ]$X4)

######################### Creation des CRHs #########################

# On récupère ici le nombre de réplicat
nb_replicas = length(unique(all_paired_structure$paire))

# nb_replicas = 2
# La boucle ira de 1 à 500
#nb_replicas = 50

for (r in 1:nb_replicas) {
  
  print(r)
  cell = r%/%2 + 1
  chr = 1
  if(r%%2==0){
    chr = 2
    cell = r/2
  }
  
  # Initialisation du premier cluster
  if(r==1){
    all_net_result <- create_bip_graphs(
      as_tibble(
        all_paired_structure%>%
          filter(
            paire == str_c(chr, sprintf("%03d", cell)),
            X4 != 1
          )%>%
          select(-c(ends_with("_c"), "paire"))%>%
          mutate(
            ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n()))
          )
      ),
      promoters_ids, 
      chr, 
      1:16, 
      3,
      sprintf("%03d", cell)
    )
  }else{
    
    structure <- as_tibble(
      all_paired_structure%>%
        filter(paire == str_c(chr, sprintf("%03d", cell)))%>%
        select(-c(ends_with("_c"), "paire"))%>%
        mutate(
          ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n()))
        )
    )
    
    print("Cluster")
    net <- create_bip_graphs(structure, promoters_ids, chr, 1:16, 3, sprintf("%03d", cell))
    
    print("Overlap")
    overlap_edge <- edge_identity_overlap(
      net, 
      all_net_result,
      0.5
    )
    
    print("Fusion")
    all_net_result <- bip_graphs_fusion(net, all_net_result, overlap_edge)
    
  }
  
}

print(all_net_result$block1$resume_fusion)

# Sauve les données dans l'objet suivant
save(all_net_result, file = "rdata/all_rda_data/all_net_result.rda")

# library(openxlsx)
# write.xlsx(
#   data.frame(
#     resume_bloc1 = all_net_result$block1$resume_fusion
#   ), 
#   file = "rdata/fusion_250_paires_block_1.xlsx"
# )


