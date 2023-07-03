

# Inspection des résultats

nb_crhs_in_clus = dim(crhs_comparation_res$sensibility_mat)[2]
nb_crh_in_struc = dim(crhs_comparation_res$sensibility_mat)[1]

nb_repete <- matrix(
  "",
  nrow = nb_crh_in_struc,
  ncol = 3,
)
for (res in seq_len(nb_crh_in_struc)){
  # sum(fireCaller_result[[res]]$FIRE_output$Mousse_cells_clus_12_indicator)
  nb_repete[res, 1] <- names(crhs_comparation_res$sensibility_mat[, 1][res])
  nb_repete[res, 2] <- names(which.max(crhs_comparation_res$sensibility_mat[res, ]))
  nb_repete[res, 3] <- max(crhs_comparation_res$sensibility_mat[res, ])
}

nb_repete = nb_repete[nb_repete[, 3] != 0, ]
plot(as.numeric(nb_repete[, 3]), type = "l")
summary(as.numeric(nb_repete[, 3]))
table(nb_repete[, 2])

fg1 = ggplot(data.frame(data = as.numeric(nb_repete[, 3])), aes(x = data)) +
  geom_boxplot()+ geom_boxplot(fill="#FC4E07")+ labs(x = "CRHs similarity index (Indice de similarité des PCis-R)")

fg2 = ggplot(data.frame(data = as.numeric(nb_repete[, 3])), aes(x=data)) +
  geom_histogram(aes(y=..density..), fill="#FC4E07")+
  geom_density(alpha=.2, fill="#E7B800")+ labs(x = "CRHs similarity index (Indice de similarité des PCis-R)", y = "Number of CRHs (Nombre de CRHs)")

figure <- ggarrange(fg1, fg2, ncol = 2, nrow = 1)
figure


s = crhs_comparation_res$sensibility_mat
summary(s)
s[s>0.5]
sum(s>0.5)
which(s>0.5)

t = apply(s, 2, function(x) length(names(which(x>0.5)))>0)
t[t==TRUE]
dim(clu_chrs_result$cluster3$crhs$crh2$mat_incidence)
dim(clu_chrs_result$cluster5$crhs$crh3$mat_incidence)

# Cluster 5 CRH1
crh = clu_chrs_result$cluster5$crhs$crh1
dim(crh$mat_incidence)
