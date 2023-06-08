# Chargement des librairies
library(tidyverse)
library(igraph)
library(ggpubr)

# Chargement des données
load("rdata/all_net_result_complex_.rda")

# Combien de contacts avons nous en moyenne dans les matrices ?

# Inspection des matrices
ncrhs = 1173

crhs_inspections = matrix(
  0,
  nrow = ncrhs,
  ncol = 4
)

ln = 1
for(bl in 2:length(all_net_result_complex_)) {
  
  for (crh in seq_len(length(all_net_result_complex_[[bl]]$crhs))) {
    if(length(all_net_result_complex_[[bl]]$crhs[[crh]])>1){
      # print(paste(bl, " ", crh, " ", length(all_net_result_complex[[bl]]$crhs[[crh]])))
      crhs_inspections[ln, 1] = nrow(all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence)
      crhs_inspections[ln, 2] = ncol(all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence)
      crhs_inspections[ln, 3] = sum(all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence)
      if(max(all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence)>1){
        print(paste(bl, " ", crh))
      }
      if((crhs_inspections[ln, 1] == 1 & crhs_inspections[ln, 2] == 1) | (crhs_inspections[ln, 1] == 2 & crhs_inspections[ln, 2] == 1) | (crhs_inspections[ln, 1] == 1 & crhs_inspections[ln, 2] == 2)){
        crhs_inspections[ln, 4] = 1
      }
      ln = ln + 1
    }
    
  }
  
}

crhs_inspections

# Combien de CRHs moins complexes ?
sum(crhs_inspections[, 4]==1)

mean(crhs_inspections[, 3])

######################################################################################"

# Inspection des nouveaux CRHs
## Initialisation du dataframe
crhs_inspections_ = data.frame(
  name =  "",
  nb_edges = NA,
  nb_vertices = NA,
  nb_sub = NA
)

for(bl in 2:length(all_net_result_complex_)) {
  
  for (crh in seq_len(length(all_net_result_complex_[[bl]]$crhs))) {
    if(length(all_net_result_complex_[[bl]]$crhs[[crh]])>1){
      
      compute_dist_bin = all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence
      net_bip <- graph_from_incidence_matrix(compute_dist_bin)
      
      # if(is_simple(net_bip)){
      #   net_bip <- simplify(net_bip, remove.multiple = FALSE, remove.loops = TRUE) 
      # }
      net_components_bip <- components(net_bip, mode = c("weak", "strong"))
      
      ln = nrow(crhs_inspections_) + 1
      if(net_components_bip$no==1){
        crhs_inspections_[ln, "name"] = paste("block", bl, "-crh", crh)
        crhs_inspections_[ln, "nb_vertices"] = length(V(net_bip))
        crhs_inspections_[ln, "nb_edges"] = length(E(net_bip))
        crhs_inspections_[ln, "nb_sub"] = net_components_bip$no
      }else{
        for (sub_crhs in seq_len(net_components_bip$no)) {
          vert_ids <- V(net_bip)[net_components_bip$membership == sub_crhs]
          sub_net <- induced_subgraph(net_bip, vert_ids)
          
          crhs_inspections_[ln, "name"] = paste("block", bl, "-crh", crh, "-sub_crh", sub_crhs)
          crhs_inspections_[ln, "nb_vertices"] = length(V(sub_net))
          crhs_inspections_[ln, "nb_edges"] = length(E(sub_net))
          crhs_inspections_[ln, "nb_sub"] = net_components_bip$no
          ln = ln + 1
        }
      }
    }
    
  }
  
}

crhs_inspections_ = crhs_inspections_[crhs_inspections_$name != "", ]
view(crhs_inspections_)

summary(crhs_inspections_$nb_edges)
summary(crhs_inspections_$nb_vertices)

crhs_inspections_tbl = data.frame(
  type = c(rep("edges", length(table(crhs_inspections_$nb_edges))), rep("vertices", length(table(crhs_inspections_$nb_vertices)))),
  n = c(data.frame(table(crhs_inspections_$nb_edges))[, 1], data.frame(table(crhs_inspections_$nb_vertices))[, 1]),
  freq = c(data.frame(table(crhs_inspections_$nb_edges))[, 2], data.frame(table(crhs_inspections_$nb_vertices))[, 2])
)

seuil_freq = 10
# width=0.7

ggplot(data=crhs_inspections_tbl[crhs_inspections_tbl$freq>seuil_freq, ], 
  aes(x=n, y=freq, fill = type)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

ggplot(data=crhs_inspections_tbl[crhs_inspections_tbl$freq>seuil_freq, ], 
       aes(x=n, y=freq, fill = type)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()

fg1 = ggplot(
  data=crhs_inspections_tbl[crhs_inspections_tbl$freq>seuil_freq & crhs_inspections_tbl$type=="edges", ], 
  aes(x=n, y=freq)) +
  geom_bar(stat="identity", position=position_dodge(), fill="#FC4E07")+
  theme_minimal()+ labs(x = "Nombre d'arrêtes", y = "Fréquences")


fg2 = ggplot(
  data=crhs_inspections_tbl[crhs_inspections_tbl$freq>seuil_freq & crhs_inspections_tbl$type=="vertices", ], 
  aes(x=n, y=freq)) +
  geom_bar(stat="identity", position=position_dodge(), fill="#E7B800")+
  theme_minimal()+ labs(x = "Nombre de noeuds", y = "Fréquences")


figure <- ggarrange(fg1, fg2, ncol = 2, nrow = 1)
figure


