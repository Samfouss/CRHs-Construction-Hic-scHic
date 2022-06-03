source("./scripts/load_save_data.R")
source("./scripts/compute_distances_fn.R")
source("./scripts/create_2D_graphics_fn.R")
source("./scripts/create_3D_graphics_fn.R")
source("./scripts/chevauche_CRHs_fn.R")
source("./scripts/find_min_by_rank.R")

######################### Creation des CRHs #########################

### Structure 1 block 1
structure_1_1_net <- create_clust_graph(structure_1, 1, 1, 3)
structure_1_1_net_comp <- components(structure_1_1_net$net, mode = c("weak", "strong"))
structure_1_1 <- data.frame(
  crh_id = structure_1_1_net_comp$membership
)%>%
  rownames_to_column('ID')%>%
  full_join(structure_1_1_net$edges, by = "ID")%>%
  full_join(structure_1[structure_1$X4==1, ], by = "ID")%>%
  relocate(X1, X2, X3, X4, ID, crh_id, nb_edges)%>%
  filter(!is.na(nb_edges))

distance_diff1 <- structure_1_1%>%
  filter(!is.na(nb_edges))%>%
  group_by(crh_id) %>%
  summarise(
    m_X1 = round(mean(X1), 2),
    m_X2 = round(mean(X2), 2),
    m_X3 = round(mean(X3), 2),
  )%>%
  relocate(starts_with("m"), crh_id)

### Structure 2 block 1
structure_2_1_net <- create_clust_graph(structure_2, 2, 1, 3)
structure_2_1_net_comp <- components(structure_2_1_net$net, mode = c("weak", "strong"))
structure_2_1 <- data.frame(
  crh_id = structure_2_1_net_comp$membership
)%>%
  rownames_to_column('ID')%>%
  full_join(structure_2_1_net$edges, by = "ID")%>%
  full_join(structure_2[structure_2$X4==1, ], by = "ID")%>%
  relocate(X1, X2, X3, X4, ID, crh_id, nb_edges)%>%
  filter(!is.na(nb_edges))

distance_diff2 <- structure_2_1%>%
  filter(!is.na(nb_edges))%>%
  group_by(crh_id) %>%
  summarise(
    m_X1 = round(mean(X1), 2),
    m_X2 = round(mean(X2), 2),
    m_X3 = round(mean(X3), 2),
  )%>%
  relocate(starts_with("m"), crh_id)

############# Représentation graphique des CRHs construis #############

### On choisi ici de representer les 4 premier dans les deux premiers replicats

grp <- plot_3D_plotly(
  structure_1_1%>%
    filter(!is.na(nb_edges))%>%
    mutate(
      rep = "1"
    )%>%
    filter(crh_id %in% 1:4)%>%
    bind_rows(
      structure_2_1%>%
        filter(!is.na(nb_edges))%>%
        mutate(
          rep = "2"
        )%>%
        filter(crh_id %in% 1:4)
    )%>%
    mutate(
      crh_id_ = paste0(rep, crh_id)
    ),
  "crh_id_"
)
grp

################### Methodes d'appariement des CRHs ###################

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

overlap_edg_max
sum(overlap_edg_max[,1]==overlap_edg_max[,2])
### Comparaison des résulats des méthodes 
dist_r1_r2_bl1_min
overlap_min
overlap_edg_max














