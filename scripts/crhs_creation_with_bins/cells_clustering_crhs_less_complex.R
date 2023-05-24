# Chargement de packages
library("tidyverse")

# construction de la matrice en retirant les CRHs non complexes
load("rdata/all_net_result.rda")
all_net_result_complex = all_net_result

for (i in 2:length(all_net_result_complex)) {
  
  for (k in seq_len(length(all_net_result_complex[[i]]$crhs))) {
    mat = all_net_result_complex[[i]]$crhs[[k]]$mat_incidence
    
    if((nrow(mat) == 1 & ncol(mat) == 1) | (nrow(mat) == 2 & ncol(mat) == 1)){
      all_net_result_complex[[i]]$crhs[[k]] <- -1
      all_net_result_complex[[i]]$resume_fusion[k] <- "-1"
    }
  }
}

########### Sauvegarde des données ###########
save(all_net_result_complex, file = "rdata/all_net_result_complex.rda")

