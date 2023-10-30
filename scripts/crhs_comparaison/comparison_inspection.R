load("rdata/all_rda_data/crhs_comparation_res.rda")

# Inspection des résultats

nb_crhs_in_clus = dim(crhs_comparation_res$specificity_mat)[2]
nb_crh_in_struc = dim(crhs_comparation_res$specificity_mat)[1]

highest_spec <- matrix(
  "",
  nrow = nb_crh_in_struc,
  ncol = 3,
)

spec_mat = crhs_comparation_res$specificity_mat
spec_mat[is.na(spec_mat)] = 0

for (res in seq_len(nb_crh_in_struc)){
  # sum(fireCaller_result[[res]]$FIRE_output$Mousse_cells_clus_12_indicator)
  highest_spec[res, 1] <- names(spec_mat[, 1][res])
  highest_spec[res, 2] <- names(which.max(spec_mat[res, ]))
  highest_spec[res, 3] <- max(spec_mat[res, ])
}

highest_spec_ = as.numeric(highest_spec[, 3])
highest_spec_ = highest_spec_[!is.na(highest_spec_)]
highest_spec_ = highest_spec_[highest_spec_ != 0]
plot(highest_spec_, type = "l")
summary(highest_spec_)

fg1 = ggplot(data.frame(data = highest_spec_), aes(x = data)) +
  geom_boxplot()+ geom_boxplot(fill="#FC4E07")+ labs(x = "Spécificité")

fg1

fg2 = ggplot(data.frame(data = highest_spec_), aes(x=data)) +
  geom_histogram(aes(y=..density..), fill="#FC4E07")+
  geom_density(alpha=.2, fill="#E7B800")+ labs(x = "Spécificité")

fg2

figure <- ggarrange(fg1, fg2, ncol = 2, nrow = 1)
figure


# Inspection des résultats

nb_crhs_in_clus = dim(crhs_comparation_res$sensibility_mat)[2]
nb_crh_in_struc = dim(crhs_comparation_res$sensibility_mat)[1]

highest_sen <- matrix(
  "",
  nrow = nb_crh_in_struc,
  ncol = 3,
)

sen_mat = crhs_comparation_res$sensibility_mat
sen_mat[is.na(sen_mat)] = 0

for (res in seq_len(nb_crh_in_struc)){
  # sum(fireCaller_result[[res]]$FIRE_output$Mousse_cells_clus_12_indicator)
  highest_sen[res, 1] <- names(sen_mat[, 1][res])
  highest_sen[res, 2] <- names(which.max(sen_mat[res, ]))
  highest_sen[res, 3] <- max(sen_mat[res, ])
}

highest_sen_ = as.numeric(highest_sen[, 3])
highest_sen_ = highest_sen_[!is.na(highest_sen_)]
highest_sen_ = highest_sen_[highest_sen_ != 0]
plot(highest_sen_, type = "l")
summary(highest_sen_)

fg1 = ggplot(data.frame(data = highest_sen_), aes(x = data)) +
  geom_boxplot()+ geom_boxplot(fill="#FC4E07")+ labs(x = "CRHs similarity index (Indice de similarité des PCis-R)")

fg2 = ggplot(data.frame(data = highest_sen_), aes(x=data)) +
  geom_histogram(aes(y=..density..), fill="#FC4E07")+
  geom_density(alpha=.2, fill="#E7B800")+ labs(x = "CRHs similarity index (Indice de similarité des PCis-R)", y = "Number of CRHs (Nombre de CRHs)")

figure <- ggarrange(fg1, fg2, ncol = 2, nrow = 1)
figure

########### Inspection des CRHs ayant une grande sensibilité au niveau des structures ##########

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

plot(
  table(as.numeric(struc_crhs_inspection[, 4])),
  xlab = "Nombre de CRHs dans les structures",
  ylab = "Nombre de répétition"
  )
plot(
  table(as.numeric(struc_crhs_inspection[, 1])),
  xlab = "Nombre de CRHs dans les structures",
  ylab = "Nombre de promoters"
)
plot(
  table(as.numeric(struc_crhs_inspection[, 2])),
  xlab = "Nombre de CRHs dans les structures",
  ylab = "Nombre d'enhancers"
)

########## Inspection des CRHs ayant une grande specificité au niveau des structures ##########

# crhs_comparation_res$specificity_mat
spec = crhs_comparation_res$specificity_mat
spec[is.na(spec)] = 0
cells_crhs = apply(spec, 1, function(x) names(which(x>0.5)))
structure_crhs = apply(spec, 2, function(x) names(which(x>0.5)))

struc_crhs_inspection = matrix(
  "",
  nrow = nrow(spec),
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

plot(
  table(as.numeric(struc_crhs_inspection[, 4])),
  xlab = "Nombre de CRHs dans les structures",
  ylab = "Nombre de répétition"
)
plot(
  table(as.numeric(struc_crhs_inspection[, 1])),
  xlab = "Nombre de CRHs dans les structures",
  ylab = "Nombre de promoters"
)
plot(
  table(as.numeric(struc_crhs_inspection[, 2])),
  xlab = "Nombre de CRHs dans les structures",
  ylab = "Nombre d'enhancers"
)





