# Exécution des scripts nécéssaires
source("./scripts/crhs_creation_with_bins/load_save_data.R")
source("./scripts/crhs_creation_with_bins/create_bip_graphs_fn.R")
source("./scripts/crhs_creation_with_bins/chevauche_CRHs_fn.R")
source("./scripts/crhs_creation_with_bins/fusionne_replicats_fn.R")

# On peut voir que dans tous blocks, à des proportions différentes, on y trouve au moins un promoter
# sort(unique(structure_1[promoters_ids, ]$X4))
# promoters_ids[which(promoters_ids<341)]
# 
# # Un petit resumé du nombre de promoters par bloc
table(structure_1[promoters_ids, ]$X4)
table(structure_2[promoters_ids, ]$X4)
table(structure_3[promoters_ids, ]$X4)

######################### Creation des CRHs #########################

# On récupère ici le nombre de réplicat
nb_replicas = length(unique(all_paired_structure$paire))

# nb_replicas = 2
# La boucle ira de 1 à 500
#nb_replicas = 50

res = "3Mb"
seuil = 0.5
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
    # print(cell)
    # print(chr)
    all_net_result <- create_bip_graphs(
      as_tibble(
        all_paired_structure%>%
          filter(paire == str_c(chr, sprintf("%03d", cell)))%>%
          select(-c(ends_with("_c"), "paire"))%>%
          mutate(
            ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n()))
          )
      ),
      promoters_ids, 
      chr, 
      1:16, 
      3,
      sprintf("%03d", cell), 
      resolution=res
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
    net <- create_bip_graphs(
      structure, 
      promoters_ids, 
      chr, 
      1:16, 
      3, 
      sprintf("%03d", cell), 
      resolution=res
    )
    
    print("Overlap")
    overlap_edge <- edge_identity_overlap(
      net, 
      all_net_result,
      seuil = seuil
    )
    
    print("Fusion")
    all_net_result <- bip_graphs_fusion(net, all_net_result, overlap_edge)
    
  }
  
}



if(res=="2Mb"){
  save(all_net_result, file = "rdata/all_rda_data/all_net_result_2Mb.rda")
}else if(res=="3Mb"){
  save(all_net_result, file = "rdata/all_rda_data/all_net_result_3Mb.rda")
}else{
  save(all_net_result, file = "rdata/all_rda_data/all_net_result_6Mb.rda")
}

load("rdata/all_rda_data/all_net_result_3Mb.rda")
ncrhs = 0
for (bl in 1:16) {
  print(length(all_net_result[[bl]]$crhs))
  ncrhs = ncrhs + length(all_net_result[[bl]]$crhs)
}
ncrhs

# print(all_net_result$block1$resume_fusion)
# library(openxlsx)
# write.xlsx(
#   data.frame(
#     resume_bloc1 = all_net_result$block1$resume_fusion
#   ), 
#   file = "rdata/fusion_250_paires_block_1.xlsx"
# )
  

