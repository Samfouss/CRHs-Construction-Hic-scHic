# Chargelent des données

load("rdata/clu_chrs_result.rda")

# Inspection des CRHS
ncrhs = 17

crhs_inspections = matrix(
  0,
  nrow = ncrhs,
  ncol = 4
)

ln = 1
for (crh in seq_len(length(clu_chrs_result))) {
  
  for (crh_ in seq_len(length(clu_chrs_result[[crh]]$crhs))) {
    
    crhs_inspections[ln, 1] = nrow(clu_chrs_result[[crh]]$crhs[[crh_]]$mat_incidence)
    crhs_inspections[ln, 2] = ncol(clu_chrs_result[[crh]]$crhs[[crh_]]$mat_incidence)
    crhs_inspections[ln, 3] = sum(clu_chrs_result[[crh]]$crhs[[crh_]]$mat_incidence)
    if((crhs_inspections[ln, 1] == 1 & crhs_inspections[ln, 2] == 1) | (crhs_inspections[ln, 1] == 2 & crhs_inspections[ln, 2] == 1) | (crhs_inspections[ln, 1] == 1 & crhs_inspections[ln, 2] == 2)){
      crhs_inspections[ln, 4] = 1
    }
    ln = ln + 1
  }
  
}

crhs_inspections

# Combien de CRHs moins complexes ?
sum(crhs_inspections[, 4]==1)



######################################################################################"

# Inspection des nouveaux CRHs
## Initialisation du dataframe
clus_crhs_inspections = data.frame(
  name =  "",
  nb_edges = NA,
  nb_vertices = NA,
  nb_sub = NA
)

for(clus in 1:length(clu_chrs_result)) {
  
  for (crh in seq_len(length(clu_chrs_result[[clus]]$crhs))) {
    compute_dist_bin = clu_chrs_result[[clus]]$crhs[[crh]]$mat_incidence
    net_bip <- graph_from_incidence_matrix(compute_dist_bin)
    
    # if(is_simple(net_bip)){
    #   net_bip <- simplify(net_bip, remove.multiple = FALSE, remove.loops = TRUE) 
    # }
    net_components_bip <- components(net_bip, mode = c("weak", "strong"))
    
    ln = nrow(clus_crhs_inspections) + 1
    if(net_components_bip$no==1){
      clus_crhs_inspections[ln, "name"] = paste("clus", clus, "-crh", crh)
      clus_crhs_inspections[ln, "nb_vertices"] = length(V(net_bip))
      clus_crhs_inspections[ln, "nb_edges"] = length(E(net_bip))
      clus_crhs_inspections[ln, "nb_sub"] = net_components_bip$no
    }else{
      for (sub_crhs in seq_len(net_components_bip$no)) {
        vert_ids <- V(net_bip)[net_components_bip$membership == sub_crhs]
        sub_net <- induced_subgraph(net_bip, vert_ids)
        
        clus_crhs_inspections[ln, "name"] = paste("clus", clus, "-crh", crh, "-sub_crh", sub_crhs)
        clus_crhs_inspections[ln, "nb_vertices"] = length(V(sub_net))
        clus_crhs_inspections[ln, "nb_edges"] = length(E(sub_net))
        clus_crhs_inspections[ln, "nb_sub"] = net_components_bip$no
        ln = ln + 1
      }
    }
  }
  
}

clus_crhs_inspections = clus_crhs_inspections[clus_crhs_inspections$name != "", ]
View(clus_crhs_inspections)

summary(clus_crhs_inspections$nb_edges)
summary(clus_crhs_inspections$nb_vertices)

clus_crhs_inspections_tbl = data.frame(
  type = c(rep("edges", length(table(clus_crhs_inspections$nb_edges))), rep("vertices", length(table(clus_crhs_inspections$nb_vertices)))),
  n = c(data.frame(table(clus_crhs_inspections$nb_edges))[, 1], data.frame(table(clus_crhs_inspections$nb_vertices))[, 1]),
  freq = c(data.frame(table(clus_crhs_inspections$nb_edges))[, 2], data.frame(table(clus_crhs_inspections$nb_vertices))[, 2])
)

seuil_freq = 0
# width=0.7

ggplot(data=clus_crhs_inspections_tbl[clus_crhs_inspections_tbl$freq>seuil_freq, ], 
       aes(x=n, y=freq, fill = type)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

ggplot(data=clus_crhs_inspections_tbl[clus_crhs_inspections_tbl$freq>seuil_freq, ], 
       aes(x=n, y=freq, fill = type)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()

fg1 = ggplot(
  data=clus_crhs_inspections_tbl[clus_crhs_inspections_tbl$freq>seuil_freq & clus_crhs_inspections_tbl$type=="edges", ], 
  aes(x=n, y=freq)) +
  geom_bar(stat="identity", position=position_dodge(), fill="#FC4E07")+
  theme_minimal()+ labs(x = "Nombre d'arrêtes", y = "Fréquences")


fg2 = ggplot(
  data=clus_crhs_inspections_tbl[clus_crhs_inspections_tbl$freq>seuil_freq & clus_crhs_inspections_tbl$type=="vertices", ], 
  aes(x=n, y=freq)) +
  geom_bar(stat="identity", position=position_dodge(), fill="#E7B800")+
  theme_minimal()+ labs(x = "Nombre de noeuds", y = "Fréquences")


figure <- ggarrange(fg1, fg2, ncol = 2, nrow = 1)
figure


load("rdata/all_rda_data/ideal_clu_chrs_result_3Mb.rda")


nb = 0
for (i in seq_len(length(clu_chrs_result))) {
  nb = nb + length(clu_chrs_result[[i]])
}

crhs_inspections = matrix(
  0,
  nrow = nb,
  ncol = 4
)

ln = 1

for(clus in 1:length(clu_chrs_result)) {
  
  for (crh in seq_len(length(clu_chrs_result[[clus]]))){
    crhs_inspections[ln, 1] = nrow(clu_chrs_result[[clus]][[crh]]$mat_incidence)
    crhs_inspections[ln, 2] = ncol(clu_chrs_result[[clus]][[crh]]$mat_incidence)
    ln = ln + 1
  }
}

summary(crhs_inspections[, 1])
summary(crhs_inspections[, 2])

load("rdata/all_rda_data/clu_chrs_result_3Mb.rda")

nb = 0
for (i in seq_len(length(clu_chrs_result))) {
  nb = nb + length(clu_chrs_result[[i]])
}

crhs_inspections = matrix(
  0,
  nrow = nb,
  ncol = 4
)

ln = 1

for(clus in 1:length(clu_chrs_result)) {
  
  for (crh in seq_len(length(clu_chrs_result[[clus]]))){
    crhs_inspections[ln, 1] = nrow(clu_chrs_result[[clus]][[crh]]$mat_incidence)
    crhs_inspections[ln, 2] = ncol(clu_chrs_result[[clus]][[crh]]$mat_incidence)
    ln = ln + 1
  }
}


summary(crhs_inspections[, 1])
summary(crhs_inspections[, 2])




