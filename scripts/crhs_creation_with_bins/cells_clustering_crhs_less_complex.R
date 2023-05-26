# Chargement de packages
library("tidyverse")

# Rapport du 22 08 2022

# construction de la matrice en retirant les CRHs non complexes
load("rdata/all_net_result.rda")
all_net_result_complex = all_net_result

for (i in 2:length(all_net_result_complex)) {
  
  for (k in seq_len(length(all_net_result_complex[[i]]$crhs))) {
    mat = all_net_result_complex[[i]]$crhs[[k]]$mat_incidence
    
    # Identification des CRHs moins complex
    if((nrow(mat) == 1 & ncol(mat) == 1) | (nrow(mat) == 2 & ncol(mat) == 1) | (nrow(mat) == 1 & ncol(mat) == 2)){
      all_net_result_complex[[i]]$crhs[[k]] <- -1
      all_net_result_complex[[i]]$resume_fusion[k] <- "-1"
    }
  }
}

# Compter le nombre de CRHs retenu

nb_crhs = 0
block_len = 0
for (bl in 2:16) {
  b = 0
  for (i in seq_len(length(all_net_result_complex[[bl]]$crhs))) {
    
    if(length(all_net_result_complex[[bl]]$crhs[[i]])>1){
      
      b = b + 1
    }
  }
  
  nb_crhs = nb_crhs + b
  block_len = block_len + b*b
}

nb_crhs
block_len

# Est ce qu'un CRHs est inclus dans l'autre après degenerescence ?
load("rdata/all_net_result_complex.rda")
# Etape 3 : fonction permettant de faire la réduction de matrices
source("scripts/crhs_comparaison/degenerationMatrix_fn.R")

crhs_struc_fusion = matrix(
  0,
  nrow = block_len,
  ncol = 4
)

l = 1

for (bl in 2:16) {
  complex_crh = vector(mode = "numeric", length = length(all_net_result_complex[[bl]]$crhs))
  index = 0
  # Identification des matrices complexes dans la liste des CRHs
  for (k in seq_len(length(all_net_result_complex[[bl]]$crhs))) {
    if(length(all_net_result_complex[[bl]]$crhs[[k]])>1){
      complex_crh[index] = k
      index = index + 1
    }
  }
  
  complex_crh = complex_crh[complex_crh != 0]
  
  for (i in complex_crh) {
    
    for (j in complex_crh) {
      
      if(i > j){
        
        # Initilisation de la premières matrice d'incidence parmi les CRHs complexes
        mat_inc_init = all_net_result_complex[[bl]]$crhs[[i]]$mat_incidence
        mat_inc_init_deg = degenerationMatrix(mat_inc_init)
        
        mat_inc = all_net_result_complex[[bl]]$crhs[[j]]$mat_incidence
        mat_inc_deg = degenerationMatrix(mat_inc)
        
        # redimenssion des deux matrices afin de calcluer les statistiques sur les intersections
        rowmat = intersect(rownames(mat_inc_deg), rownames(mat_inc_init_deg))
        colmat = intersect(colnames(mat_inc_deg), colnames(mat_inc_init_deg))
        
        mat_deg_redim1 <- matrix(
          0,
          ncol=length(colmat), 
          nrow=length(rowmat), 
          dimnames=list(rowmat, colmat)
        )
        mat_deg_redim2 <- mat_deg_redim1
        
        mat_deg_redim1 <- mat_inc_deg[rowmat, colmat]
        mat_deg_redim2 <- mat_inc_init_deg[rowmat, colmat]
        
        if((min(sum(mat_deg_redim1==1), sum(mat_deg_redim2==1))) != 0){
          intersection = sum(mat_deg_redim1[mat_deg_redim2==mat_deg_redim1]==1)/(min(sum(mat_deg_redim1==1), sum(mat_deg_redim2==1)))
          # Sauvegarde des données
          crhs_struc_fusion[l, 1] = bl
          crhs_struc_fusion[l, 2] = i
          crhs_struc_fusion[l, 3] = j
          crhs_struc_fusion[l, 4] = intersection
          l = l + 1
        }
      }
      
    }
  }
  
}

crhs_struc_fusion = crhs_struc_fusion[crhs_struc_fusion[, 1]!= 0, ]

summary(crhs_struc_fusion[, 4])

save(crhs_struc_fusion, file = "rdata/crhs_struc_fusion.rda")
########### Sauvegarde des données ###########
save(all_net_result_complex, file = "rdata/all_net_result_complex.rda")

