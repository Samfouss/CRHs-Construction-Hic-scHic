

#' Cette fonction donne une representation graphique des points en 2D en utilisant
#' ggplot
#'
#' @description
#' `plot_3D_plotly` Cette fonction donne une representation graphique des 
#' points en 2D en utilisant ggplot
#'
#' @param data_to_plot 
#' 

plot_2D <- function(data_too_plot){
  
  g1 <- ggplot(data_too_plot, aes(X1, X2))+
    geom_point()+
    geom_point(aes(colour = factor(replicat)))
  
  g2 <- ggplot(data_too_plot, aes(X1, X3))+
    geom_point()+
    geom_point(aes(colour = factor(replicat)))
  
  g3 <- ggplot(data_too_plot, aes(X2, X3))+
    geom_point()+
    geom_point(aes(colour = factor(replicat)))
  
  ggarrange(g1, g2, g3 + rremove("x.text"), 
            labels = c("A", "B", "C"),
            ncol = 2, nrow = 2)
}


