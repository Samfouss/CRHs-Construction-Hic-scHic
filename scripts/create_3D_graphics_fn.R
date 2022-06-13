
#' Cette fonction donne une representation graphique des points en 3D en utilisant
#' plotly
#'
#' @description
#' `plot_3D_plotly` Cette fonction donne une representation graphique des 
#' points en 3D en utilisant plotly
#'
#' @param data_rep1_to_plot : le premier replicat à representer
#' @param data_rep1_to_plot : le deuxième replicat à representer
#' @param block_to_rep : le numero du block à representer
#' @param crh_match_list : la liste contenant l'ensemble des CRHs qui matchent pour tous les blocks

# sprintf("%03d", block_to_rep)
plot_3D_plotly <- function(data_rep1_to_plot, data_rep2_to_plot, block_to_rep, crh_match_list){
  
  crh_match_list1 = crh_match_list[[block_to_rep]]$crh_edge_match[1, ]
  crh_match_list2 = crh_match_list[[block_to_rep]]$crh_edge_match[2, ]
  
  data_to_plot <- data_rep1_to_plot%>%
    filter(X4==block_to_rep)%>%
    mutate(
      rep_num = "1"
    )%>%
    filter(crh_id %in% crh_match_list1)%>%
    bind_rows(
      data_rep2_to_plot%>%
        filter(X4==block_to_rep)%>%
        mutate(
          rep_num = "2"
        )%>%
        filter(crh_id %in% crh_match_list2)
    )%>%
    mutate(
      crh_id_ = paste0(rep_num, sprintf("%02d", crh_id))
    )
  
  data_to_plot
  
  data_to_plot <- data_to_plot %>%
    group_by(crh_id_) %>%
    mutate(group_id =cur_group_id())%>%
    ungroup()
  
  # distinctColorPalette(length(unique(data_to_plot$group_id)))
  data_to_plot <- data_to_plot%>%
    left_join(
      data.frame(color = distinctColorPalette(length(unique(data_to_plot$group_id))))%>%
        rownames_to_column("group_id")%>%
        mutate(group_id = as.numeric(group_id)),
      by = "group_id"
    )%>%
    arrange(group_id)%>%
    mutate_(
      group_id = "crh_id_"
    )
  
  fig <- data_to_plot%>%
    plot_ly(
      x = ~X1, 
      y = ~X2, 
      z = ~X3, 
      color = ~group_id, 
      colors = unique(unlist(data_to_plot[, "color"]))
    )
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(
    scene = list(xaxis = list(title = 'X1'),
                 yaxis = list(title = 'X2'),
                 zaxis = list(title = 'X3'))
  )
  fig
  
}


