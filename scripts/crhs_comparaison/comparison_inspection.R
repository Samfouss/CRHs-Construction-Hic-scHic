

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

############################# Inspection des CRHs ayant une grande sensibilité au niveau des structures

# crhs_comparation_res$specificity_mat
sen = crhs_comparation_res$sensibility_mat
sen[is.na(sen)] = 0
cells_crhs = apply(sen, 1, function(x) names(which(x>0.5)))
structure_crhs = apply(sen, 2, function(x) names(which(x>0.5)))

struc_crhs_inspection = matrix(
  "",
  nrow = nrow(sen),
  ncol = 4
)

l = 0
for (i in seq_len(length(structure_crhs))) {
  if(length(structure_crhs[[i]])>0){
    for (j in seq_len(length(structure_crhs[[i]]))) {
      bl = as.numeric(unlist(strsplit(structure_crhs[[i]][[j]][1], "[_]"))[2])
      crh = as.numeric(unlist(strsplit(structure_crhs[[i]][[j]][1], "[_]"))[4])
      mat = all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence
      crh_name = paste0("block ", sprintf("%02d", bl), "- CRH ", crh)
      
      if(crh_name %in% struc_crhs_inspection[, 3]){
        r = which(struc_crhs_inspection[, 3]==crh_name)
        rep = as.numeric(struc_crhs_inspection[r, 4]) + 1
        struc_crhs_inspection[r, 4] =  as.character(rep)
      }else{
        struc_crhs_inspection[l, 1] =  as.character(dim(mat)[1])
        struc_crhs_inspection[l, 2] =  as.character(dim(mat)[2])
        struc_crhs_inspection[l, 3] =  crh_name
        struc_crhs_inspection[l, 4] = as.character(1)
        l = l + 1
      }
    }
  }
}

struc_crhs_inspection = struc_crhs_inspection[struc_crhs_inspection[, 3]!="", ]

save(
  struc_crhs_inspection, 
  file = "rdata/all_rda_data/struc_crhs_inspection.rda"
)



