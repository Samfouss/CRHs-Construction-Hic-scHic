library(ggplot2)

load("rdata/all_rda_data/crhs_comparaison_res.rda")

load("rdata/all_rda_data/ideal_crhs_comparaison_res.rda")
# Inspection des résultats

nb_crhs_in_clus = dim(crhs_comparaison_res$specificity_mat)[2]
nb_crh_in_struc = dim(crhs_comparaison_res$specificity_mat)[1]

highest_spec <- matrix(
  "",
  nrow = nb_crh_in_struc,
  ncol = 3,
)

spec_mat = crhs_comparaison_res$specificity_mat

l = 1
for (res in seq_len(nb_crh_in_struc)){
  # sum(fireCaller_result[[res]]$FIRE_output$Mousse_cells_clus_12_indicator)
  if(!all(is.na(spec_mat[res, ]))){
    highest_spec[l, 1] <- names(spec_mat[, 1][res])
    highest_spec[l, 2] <- names(which.max(spec_mat[res, ]))
    highest_spec[l, 3] <- max(spec_mat[res, ][!is.na(spec_mat[res, ])])
    l = l + 1
  }
}

highest_spec = highest_spec[highest_spec[, 1] != "", ]

highest_spec_ = as.numeric(highest_spec[, 3])
# highest_spec_ = highest_spec_[!is.na(highest_spec_)]
# highest_spec_ = highest_spec_[highest_spec_ != 0]
length(highest_spec_)
plot(
  highest_spec_, 
  type = "l",
  ylab = "Spécificité",
  xlab = "Nombre de CRHs"
)
hist(
  highest_spec_,
  xlab = "Spécificité",
  ylab = "Nombre de CRHs",
  main = ""
)
boxplot(highest_spec_, horizontal = TRUE)
summary(highest_spec_)

# ggplot(data.frame(data = highest_spec_), aes(x = data)) +
#   geom_boxplot()+ geom_boxplot(fill="#FC4E07")+ labs(x = "Spécificité", y = "Number of CRHs (Nombre de CRHs)")
# 
# ggplot(data.frame(data = highest_spec_), aes(x=data)) +
#   geom_histogram(aes(y=..density..), fill="#FC4E07")+
#   geom_density(alpha=.2, fill="#E7B800")+ labs(x = "Spécificité", y = "Number of CRHs (Nombre de CRHs)")


# Inspection des résultats

nb_crhs_in_clus = dim(crhs_comparaison_res$sensibility_mat)[2]
nb_crh_in_struc = dim(crhs_comparaison_res$sensibility_mat)[1]

highest_sen <- matrix(
  "",
  nrow = nb_crh_in_struc,
  ncol = 3,
)

sen_mat = crhs_comparaison_res$sensibility_mat

l = 1
for (res in seq_len(nb_crh_in_struc)){
  # sum(fireCaller_result[[res]]$FIRE_output$Mousse_cells_clus_12_indicator)
  if(!all(is.na(sen_mat[res, ]))){
    highest_sen[l, 1] <- names(sen_mat[, 1][res])
    highest_sen[l, 2] <- names(which.max(sen_mat[res, ]))
    highest_sen[l, 3] <- max(sen_mat[res, ][!is.na(sen_mat[res, ])])
    l = l + 1
  }
}

highest_sen = highest_sen[highest_sen[, 1] != "", ]

highest_sen_ = as.numeric(highest_sen[, 3])
# highest_sen_ = highest_sen_[!is.na(highest_sen_)]
length(highest_sen_)
# highest_sen_ = highest_sen_[highest_sen_ != 0]
plot(
  highest_sen_, 
  type = "l",
  ylab = "Sensibilité",
  xlab = "Nombre de CRHs",
  #col="gray"
)
hist(
  highest_sen_,
  xlab = "Sensibilité",
  ylab = "Nombre de CRHs",
  main = "",
  #col="gray"
)
boxplot(highest_sen_, horizontal = TRUE)
summary(highest_sen_)

# Calcul de li'intervalle de confiance
## Borne inférieure

# ggplot(data.frame(data = highest_sen_), aes(x = data)) +
#   geom_boxplot()+ geom_boxplot(fill="#FC4E07")+ labs(x = "sensibilité", y = "Number of CRHs (Nombre de CRHs)")
# 
# ggplot(data.frame(data = highest_sen_), aes(x=data)) +
#   geom_histogram(aes(y=..density..), fill="#FC4E07")+
#   geom_density(alpha=.2, fill="#E7B800")+ labs(x = "sensibilité", y = "Number of CRHs (Nombre de CRHs)")

################## Au niveau des cellules ##################

########### Inspection des CRHs ayant une grande sensibilité au niveau des structures ##########

sen = crhs_comparaison_res$sensibility_mat
sen[is.na(sen)] = 0
cells_crhs = apply(sen, 1, function(x) names(which(x>=0.5)))

clusters_name = c()
for (i in seq_len(length(cells_crhs))) {
  if(length(cells_crhs[[i]])>0){
    clusters_name = c(clusters_name, cells_crhs[[i]])
  }
}

hist(table(clusters_name), xlab = "Nombre de répetition des CRHs dans les cluster", main = "Les CRHs dans les clusters dont la sensibilié est superieure ou égale à 0.5")

cells_crhs_inspection = matrix(
  "",
  nrow = ncol(sen),
  ncol = 4
)

l = 0
for (i in seq_len(length(cells_crhs))) {
  if(length(cells_crhs[[i]])>0){
    for (j in seq_len(length(cells_crhs[[i]]))) {
      clus = as.numeric(unlist(strsplit(cells_crhs[[i]][[j]][1], "[_]"))[2])
      
      crh = as.numeric(unlist(strsplit(cells_crhs[[i]][[j]][1], "[_]"))[4])
      clu_chrs_result[[clus]][[j]]$crhs[[i]]$mat_incidence
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
spec = crhs_comparaison_res$specificity_mat
spec[is.na(spec)] = 0
cells_crhs = apply(spec, 1, function(x) names(which(x>0.5)))

clusters_name = c()
for (i in seq_len(length(cells_crhs))) {
  if(length(cells_crhs[[i]])>0){
    clusters_name = c(clusters_name, cells_crhs[[i]])
  }
}

hist(table(clusters_name), xlab = "Nombre de répetition des CRHs dans les cluster", main = "Les CRHs dans les clusters dont la specificité est superieure ou égale à 0.5")

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

all(crhs_comparaison_res$sensibility_mat==crhs_comparaison_res2$sensibility_mat)
