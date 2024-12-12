############# Ajout des statistiques sur la complexité des réseaux
load("rdata/all_rda_data/nb_prom_enhan_matrix.rda")
# On retire les NA dans les données
nb_prom_enhan_matrix <- nb_prom_enhan_matrix[!is.na(nb_prom_enhan_matrix[, 1]), ]

nb_prom_enhan_matrix <- as.data.frame(nb_prom_enhan_matrix)
colnames(nb_prom_enhan_matrix) <- c("fold", "Promoteurs", "Enhancers")
summary(nb_prom_enhan_matrix$Promoteurs)
summary(nb_prom_enhan_matrix$Enhancers)

summary(nb_prom_enhan_matrix[nb_prom_enhan_matrix$fold==1, ]$Promoteurs)
summary(nb_prom_enhan_matrix[nb_prom_enhan_matrix$fold==1, ]$Enhancers)

summary(nb_prom_enhan_matrix[nb_prom_enhan_matrix$fold==2, ]$Promoteurs)
summary(nb_prom_enhan_matrix[nb_prom_enhan_matrix$fold==2, ]$Enhancers)

summary(nb_prom_enhan_matrix[nb_prom_enhan_matrix$fold==3, ]$Promoteurs)
summary(nb_prom_enhan_matrix[nb_prom_enhan_matrix$fold==3, ]$Enhancers)

summary(nb_prom_enhan_matrix[nb_prom_enhan_matrix$fold==4, ]$Promoteurs)
summary(nb_prom_enhan_matrix[nb_prom_enhan_matrix$fold==4, ]$Enhancers)

summary(nb_prom_enhan_matrix[nb_prom_enhan_matrix$fold==5, ]$Promoteurs)
summary(nb_prom_enhan_matrix[nb_prom_enhan_matrix$fold==5, ]$Enhancers)

