# Chargement de packages
library("FactoMineR")
library("factoextra")

cells_data <- t(cells_matrix)
cells_data[cells_data == 1] <- "in"
cells_data[cells_data == 0] <- "Not in"

cells_data <- data.frame(cells_data)

res.mca <- MCA(cells_data, graph = FALSE)
print(res.mca)
eig.val <- get_eigenvalue(res.mca)

fviz_mca_biplot(res.mca, repel = TRUE, ggtheme = theme_minimal())
fviz_mca_var(res.mca, choice = "mca.cor", repel = TRUE)
fviz_mca_var(res.mca, repel = TRUE, ggtheme= theme_minimal())
fviz_mca_var(res.mca, choice = "quanti.sup", ggtheme = theme_minimal())
fviz_mca_ind(res.mca, label = "ind.sup", ggtheme = theme_minimal())

