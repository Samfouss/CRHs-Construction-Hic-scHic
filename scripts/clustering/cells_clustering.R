# Chargement de packages
library("tidyverse")
library("FactoMineR")
library("factoextra")
library("ggpubr")

load("rdata/cells_matrix.rda")
cells_data <- cells_matrix[[1]]
for (b in 2:16) {
  cells_data <- cbind(
    cells_data,
    cells_matrix[[b]] 
  )
}

#cells_data <- t(cells_matrix)
cells_data[cells_data == 1] <- "in"
cells_data[cells_data == 0] <- "Not in"

#cells_data <- data.frame(cells_data)

res.mca <- MCA(cells_data, graph = FALSE)
print(res.mca)
eig.val <- get_eigenvalue(res.mca)

fviz_mca_biplot(res.mca, repel = TRUE, ggtheme = theme_minimal())
fviz_mca_var(res.mca, choice = "mca.cor", repel = TRUE)
fviz_mca_var(res.mca, repel = TRUE, ggtheme= theme_minimal())
fviz_mca_var(res.mca, choice = "quanti.sup", ggtheme = theme_minimal())
fviz_mca_ind(res.mca, label = "ind.sup", ggtheme = theme_minimal())
fviz_mca_ind(res.mca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, ggtheme = theme_minimal())
hcpc_cluster <- HCPC(res.mca, kk=Inf, min = 3, graph = FALSE)
# Le dendogram
fviz_dend(hcpc_cluster, cex = 0.7, palette = "jco", rect = TRUE, rect_fill = TRUE, rect_border = "jco")
# Graphique
fviz_cluster(hcpc_cluster, repel = TRUE, geom = "point", main = "Classification des cellules")
table(hcpc_cluster$data.clust$clust)
hcpc_cluster$call$t


