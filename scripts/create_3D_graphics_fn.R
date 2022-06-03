
#' Cette fonction donne une representation graphique des points en 3D en utilisant
#' plotly
#'
#' @description
#' `plot_3D_plotly` Cette fonction donne une representation graphique des 
#' points en 3D en utilisant plotly
#'
#' @param data_to_plot 
#' @param color_ref 

plot_3D_plotly <- function(data_to_plot, color_ref){
  
  data_to_plot <- data_to_plot %>% 
    group_by_(color_ref) %>% 
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
      group_id = color_ref
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


#' Cette fonction donne une representation graphique des points en 3D en utilisant
#' plot3d
#'
#' @description
#' `plot_3D` Cette fonction donne une representation graphique des 
#' points en 3D en utilisant plot3d
#'
#' @param data_to_plot 
#' @param color_ref 

plot_3D <- function(data_too_plot, color_ref , display = FALSE, radus = .3, axe_rotate = c(1,0,0), rotate_speed = 3, duration = 10, gif=FALSE){
  
  data_too_plot <- data_too_plot %>% 
    group_by_(color_ref) %>% 
    mutate(group_id =cur_group_id())%>%
    ungroup()
  
  data_too_plot <- data_too_plot%>%
    left_join(
      data.frame(color = distinctColorPalette(length(unique(data_too_plot$group_id))))%>%
        rownames_to_column("group_id")%>%
        mutate(group_id = as.numeric(group_id)),
      by = "group_id"
    )
  
  par3d(windowRect = c(0, 0, 1100, 700))
  plot3d(
    x = as.numeric(unlist(data_too_plot[,1])), 
    y = as.numeric(unlist(data_too_plot[,2])), 
    z = as.numeric(unlist(data_too_plot[,3])), 
    col = data_too_plot$color, 
    type = "s", 
    radius = radus,
    xlab = "X", 
    ylab = "Y", 
    zlab = "Z"
  )
  
  if (display){
    play3d(spin3d(axis = axe_rotate, rpm = rotate_speed), duration = duration) 
  }
  
  if(gif){
    play3d(spin3d(axis = axe_rotate, rpm = rotate_speed), duration = duration)
    movie3d(
      spin3d(axis = axe_rotate, rpm = rotate_speed),
      duration = duration,
      type = "gif",
      dir = ".",
      fps = 20
    )
  }
}