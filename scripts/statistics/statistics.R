library(tidyverse)
library(ggcharts)
library(boot)

load("rdata/all_rda_data/all_net_result_complex_3Mb_.rda")
load("rdata/all_rda_data/clu_chrs_result_3Mb.rda")
load("rdata/all_rda_data/cluster_result.rda")

blocs_inspection_ = matrix(
  0,
  nrow = 16,
  ncol = 3
)
blocs_inspection__ = matrix(
  0,
  nrow = 16,
  ncol = 3
)

nb_bins = c(342, 133, 207, 179, 101, 115, 174,  93, 123, 209, 90,  99, 115,  84,  90,  94)

for (bl in 1:16) {
  b = 0
  for (i in seq_len(length(all_net_result_complex_[[bl]]$crhs))) {
    if(sum(all_net_result_complex_[[bl]]$crhs[[i]]$mat_incidence) != -1){
      b = b + 1
    }
  }
  blocs_inspection_[bl, 1] = bl
  blocs_inspection_[bl, 2] = b
  blocs_inspection_[bl, 3] = 1
  
  blocs_inspection__[bl, 1] = bl
  blocs_inspection__[bl, 2] = nb_bins[bl]
  blocs_inspection__[bl, 3] = 2
}

blocs_inspection = rbind(blocs_inspection_[blocs_inspection_[, 1]!=1, ], blocs_inspection__[blocs_inspection__[, 1]!=1, ])

blocs_inspection = tibble(
  Blocs = blocs_inspection[, 1],
  Nombre = blocs_inspection[, 2],
  type = blocs_inspection[, 3],
)%>%
  mutate(
    label = case_when(
      type == 1  ~ "Nombre de perles",
      type == 2  ~ "Nombre de réseaux"
    )
  )

ggplot(data=blocs_inspection, aes(x=Blocs, y=Nombre, fill=label)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=Nombre), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()





clusters = tibble(
  groupe = as.data.frame(table(cells_clusters$cluster))[, 1],
  taille = as.data.frame(table(cells_clusters$cluster))[, 2]
)%>%
  mutate(
    taille_ = case_when(
      taille < 40  ~ "Moins de 40",
      taille < 80  ~ "Entre 40 et 80",
      taille < 120 ~ "Entre 80 et 120",
      taille < 160 ~ "Entre 120 et 160",
      taille <= 200 ~ "Plus de 160",
    ),
    .keep = "used"
  )

ggplot(data = clusters) + 
  geom_bar(mapping = aes(x = taille_)) +
  labs(
    x = "La taille des groupes de cellules",
    y = "Nombre de groupes de cellules"
  )

clusters%>%
  summarise(
    x = 1,
    min = min(taille),
    first_Qu = quantile(taille, probs = 0.25),
    median = median(taille),
    third_Qu = quantile(taille, probs = 0.75),
    max = max(taille)
  )%>%
  ggplot(aes(x)) +
  geom_boxplot(
    aes(ymin = min, lower = first_Qu, middle = median, upper = third_Qu, ymax = max),
    stat = "identity"
  ) +
  labs(
    x = "",
    y = "La distribution de la taille des groupes de cellules"
  )

summary(clusters$taille)

######### Calcul des intervalles de confiance cas idéal ###########

load("rdata/all_rda_data/ideal_crhs_comparaison_res.rda")


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

n = length(highest_spec_)
p = median(highest_spec_)
spec_IC_inf = p - 1.96*sqrt(p*(1-p)/n)
spec_IC_inf
spec_IC_sup = p + 1.96*sqrt(p*(1-p)/n)
spec_IC_sup


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
    highest_sen[res, 1] <- names(sen_mat[, 1][res])
    highest_sen[res, 2] <- names(which.max(sen_mat[res, ]))
    highest_sen[l, 3] <- max(sen_mat[res, ][!is.na(sen_mat[res, ])])
    l = l + 1
  }
}

highest_sen = highest_sen[highest_sen[, 1] != "", ]
highest_sen_ = as.numeric(highest_sen[, 3])

summary(highest_sen_)
summary(highest_spec_)

dim(crhs_comparaison_res$specificity_mat)
dim(crhs_comparaison_res$sensibility_mat)

n = length(highest_sen_)
p = median(highest_sen_)
spec_IC_inf = p - 1.96*sqrt(p*(1-p)/n)
spec_IC_inf
spec_IC_sup = p + 1.96*sqrt(p*(1-p)/n)
spec_IC_sup

######### Calcul des intervalles de confiance avec 250 groupes ###########

load("rdata/all_rda_data/crhs_comparaison_res.rda")

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

n = length(highest_spec_)
p = median(highest_spec_)
spec_IC_inf = p - 1.96*sqrt(p*(1-p)/n)
spec_IC_inf
spec_IC_sup = p + 1.96*sqrt(p*(1-p)/n)
spec_IC_sup


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
    highest_sen[res, 1] <- names(sen_mat[, 1][res])
    highest_sen[res, 2] <- names(which.max(sen_mat[res, ]))
    highest_sen[l, 3] <- max(sen_mat[res, ][!is.na(sen_mat[res, ])])
    l = l + 1
  }
}

highest_sen = highest_sen[highest_sen[, 1] != "", ]
highest_sen_ = as.numeric(highest_sen[, 3])

n = length(highest_sen_)
p = median(highest_sen_)
spec_IC_inf = p - 1.96*sqrt(p*(1-p)/n)
spec_IC_inf
spec_IC_sup = p + 1.96*sqrt(p*(1-p)/n)
spec_IC_sup

############################################################# Résultats de la validation croisée

library(boot)
load("rdata/all_rda_data/res_matrix.rda")

all_res = res_matrix[row.names(res_matrix)=="" & !is.na(res_matrix[, 1]), ]
all_res

summary(all_res[, 3])
summary(all_res[, 4])



boot.ci(boot_res, type=c("norm", "basic", "perc"), index = 3)
boot.ci(boot_res, type=c("norm", "basic", "perc"), index = 4)


mean(boot_res$t[, 3])
mean(boot_res$t[, 4])



load("rdata/all_rda_data/ideal_boot_res.rda")

mean(boot_res$t[, 3])
mean(boot_res$t[, 4])

2*boot_res$t0 - colMeans(boot_res$t)
boot.ci(boot_res, type=c("norm", "basic", "perc"), index = 3)
boot.ci(boot_res, type=c("norm", "basic", "perc"), index = 4)








