source("./scripts/load_save_data.R")
source("./scripts/compute_distances_fn.R")
source("./scripts/create_2D_graphics_fn.R")
source("./scripts/create_3D_graphics_fn.R")
source("./scripts/chevauche_CRHs_fn.R")

####### Creation des CRHs
structure_1_1_net <- create_clust_graph(structure_1, 1, 1, 3)
structure_1_1_net_comp <- components(structure_1_1_net$net, mode = c("weak", "strong"))
# Il y a 8 qui sont des singletons
structure_1_1 <- data.frame(
  crh_id = structure_1_1_net_comp$membership
)%>%
  rownames_to_column('ID')%>%
  full_join(structure_1_1_net$edges, by = "ID")%>%
  full_join(structure_1[structure_1$X4==1, ], by = "ID")%>%
  relocate(X1, X2, X3, X4, ID, crh_id, nb_edges)

distance_diff1 <- structure_1_1%>%
  group_by(crh_id) %>%
  summarise(
    m_X1 = round(mean(X1), 2),
    m_X2 = round(mean(X2), 2),
    m_X3 = round(mean(X3), 2),
  )%>%
  relocate(starts_with("m"), crh_id)

structure_2_1_net <- create_clust_graph(structure_2, 2, 1, 3)
structure_2_1_net_comp <- components(structure_2_1_net$net, mode = c("weak", "strong"))
structure_2_1 <- data.frame(
  crh_id = structure_2_1_net_comp$membership
)%>%
  rownames_to_column('ID')%>%
  full_join(structure_2_1_net$edges, by = "ID")%>%
  full_join(structure_2[structure_2$X4==1, ], by = "ID")%>%
  relocate(X1, X2, X3, X4, ID, crh_id, nb_edges)

distance_diff2 <- structure_2_1%>%
  group_by(crh_id) %>%
  summarise(
    m_X1 = round(mean(X1), 2),
    m_X2 = round(mean(X2), 2),
    m_X3 = round(mean(X3), 2),
  )%>%
  relocate(starts_with("m"), crh_id)


###### Représentation graphique #####
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

######## Calcul des distances entre les CRHs #######
dist_r1_r2_bl1 <- m_dist_func(distance_diff1, distance_diff2)

dist_r1_r2_bl1_min <- as.matrix(
  apply(dist_r1_r2_bl1, MARGIN = 1, FUN = min)
)

dist_r1_r2_bl1_min







