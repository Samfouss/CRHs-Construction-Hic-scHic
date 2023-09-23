# Chargement de packages
library("tidyverse")

# Rapport du 22 08 2022

# construction de la matrice en retirant les CRHs non complexes
load("rdata/all_rda_data/all_net_result.rda")
all_net_result_complex = all_net_result

for (i in 2:length(all_net_result_complex)) {
  for (k in seq_len(length(all_net_result_complex[[i]]$crhs))) {
    mat = all_net_result_complex[[i]]$crhs[[k]]$mat_incidence
    
    # Identification des CRHs moins complex
    if((nrow(mat) == 1 & ncol(mat) == 1) | (nrow(mat) == 2 & ncol(mat) == 1) | (nrow(mat) == 1 & ncol(mat) == 2)){
      all_net_result_complex[[i]]$crhs[[k]]$mat_incidence <- -1
      all_net_result_complex[[i]]$more_info <- "-1"
    }
  }
}

# Compter le nombre de CRHs retenu
nb_crhs = 0
block_len = 0
for (bl in 2:16) {
  b = 0
  for (i in seq_len(length(all_net_result_complex[[bl]]$crhs))) {
    if(sum(all_net_result_complex[[bl]]$crhs[[i]]$mat_incidence) != -1){
      b = b + 1
    }
  }
  print(b)
  nb_crhs = nb_crhs + b
  block_len = block_len + b*b
}

nb_crhs
block_len
save(all_net_result_complex, file ="rdata/all_rda_data/all_net_result_complex.rda")

# Est ce qu'un CRHs est inclus dans l'autre après degenerescence ?
load("rdata/all_rda_data/all_net_result_complex.rda")

# Etape 3 : fonction permettant de faire la réduction de matrices
source("scripts/crhs_comparaison/degenerationMatrix_fn.R")

crhs_struc_fusion = matrix(
  0,
  nrow = block_len,
  ncol = 5
)

l = 1

for (bl in 2:16) {
  complex_crh = vector(mode = "numeric", length = length(all_net_result_complex[[bl]]$crhs))
  index = 0
  # Identification des matrices complexes dans la liste des CRHs
  for (k in seq_len(length(all_net_result_complex[[bl]]$crhs))) {
    if(sum(all_net_result_complex[[bl]]$crhs[[k]]$mat_incidence) != -1){
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
        rowmat = union(rownames(mat_inc_deg), rownames(mat_inc_init_deg))
        colmat = union(colnames(mat_inc_deg), colnames(mat_inc_init_deg))
        
        mat_deg_redim1 <- matrix(
          0,
          ncol=length(colmat), 
          nrow=length(rowmat), 
          dimnames=list(rowmat, colmat)
        )
        mat_deg_redim2 <- mat_deg_redim1
        
        indxA <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat_inc_deg), colnames(mat_inc_deg), FUN=paste)
        indxB <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat_inc_init_deg), colnames(mat_inc_init_deg), FUN=paste)
        mat_deg_redim1[indxA] <- mat_inc_deg
        mat_deg_redim2[indxB] <- mat_inc_init_deg
        
        if((min(sum(mat_deg_redim1==1), sum(mat_deg_redim2==1))) != 0){
          intersection = sum(mat_deg_redim1[mat_deg_redim2==mat_deg_redim1]==1)/(min(sum(mat_deg_redim1==1), sum(mat_deg_redim2==1)))
          # Sauvegarde des données
          crhs_struc_fusion[l, 1] = bl
          crhs_struc_fusion[l, 2] = i
          crhs_struc_fusion[l, 3] = j
          crhs_struc_fusion[l, 4] = intersection
          crhs_struc_fusion[l, 5] = i
          if(min(sum(mat_deg_redim1==1), sum(mat_deg_redim2==1))==sum(mat_deg_redim1==1)){
            crhs_struc_fusion[l, 5] = j
          }
          l = l + 1
        }
      }
      
    }
  }
  
}

crhs_struc_fusion = crhs_struc_fusion[crhs_struc_fusion[, 1]!= 0, ]

summary(crhs_struc_fusion[, 4])

########### Sauvegarde des résulats sur les paires de CRHs qui se chevauchent ##########
save(crhs_struc_fusion, file = "rdata/all_rda_data/crhs_struc_fusion.rda")


################################ Fusion des CRHs qui se chevauchent ##########################################"


########### Sauvegarde des résulats sur les paires de CRHs qui se chevauchent ##########
load("rdata/all_rda_data/crhs_struc_fusion.rda")

all_net_result_complex_ = all_net_result_complex
unique_crhs = unique(crhs_struc_fusion[, 1])

# Initilisation des données sur les CRHs
for (bl in unique_crhs) {
  block = crhs_struc_fusion[crhs_struc_fusion[, 1]==bl, ]
  
  crh_to_remove = sort(unique(block[block[, 4]==1, 5]))
  for (i in crh_to_remove) {
    all_net_result_complex_[[bl]]$crhs[[i]]$mat_incidence = -1
  }
}


nb_crhs = 0
for (bl in 2:16) {
  b = 0
  for (i in seq_len(length(all_net_result_complex_[[bl]]$crhs))) {
    
    if(sum(all_net_result_complex_[[bl]]$crhs[[i]]$mat_incidence) != -1){
      all_net_result_complex_[[bl]]$crhs[[i]]$mat_incidence = degenerationMatrix(all_net_result_complex_[[bl]]$crhs[[i]]$mat_incidence)
      b = b + 1
    }
  }
  # print(b)
  nb_crhs = nb_crhs + b
}
nb_crhs

# for (bl in 2:length(all_net_result_complex)) {
#   
#   # Ici on recupere
#   block = crhs_struc_fusion[crhs_struc_fusion[, 1]==bl, ]
#   unique_crhs = unique(block[, 2])  
#   for (crh in unique_crhs) {
#     
#     # On recupere ici les informations du CRHs courant à fusioner
#     all_net_result_complex_[[bl]]$crhs[[crh]] = all_net_result_complex[[bl]]$crhs[[crh]]
#     # On prend ici les informations sur le CRH fixé
#     crh_fixed = block[block[, 2]==crh, , drop=FALSE]
#     # On prend ici les ou le CRH(s) qui chevauche le mieux le CRH fixé
#     max_intersect = which(crh_fixed[, 4]==1)
#     
#     if(length(max_intersect)>0){
#       crhs_to_merge = crh_fixed[max_intersect, 3]
#       
#       for (crh_to_merge in crhs_to_merge) {
#         # if(){
#         #   
#         # }
#         # print(paste(crh_to_merge, " ", bl))
#         mat_init = all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence
#         mat_to_merge = all_net_result_complex[[bl]]$crhs[[crh_to_merge]]$mat_incidence
#         
#         rowmat = union(rownames(mat_init), rownames(mat_to_merge))
#         colmat = union(colnames(mat_init), colnames(mat_to_merge))
#         
#         mat_init_redim <- matrix(
#           0,
#           ncol=length(colmat), 
#           nrow=length(rowmat), 
#           dimnames=list(rowmat, colmat)
#         )
#         mat_to_merge_redim <- mat_init_redim
#         
#         indxA <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat_init), colnames(mat_init), FUN=paste)
#         indxB <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat_to_merge), colnames(mat_to_merge), FUN=paste)
#         mat_init_redim[indxA] <- mat_init
#         mat_to_merge_redim[indxB] <- mat_to_merge
#         
#         mat_fusion = mat_init_redim + mat_to_merge_redim
#         mat_fusion[mat_fusion>1] <- 1
#         all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence <- mat_fusion
#       } 
#     }
#     
#   }
#   
# }

nb_crhs = 0
for (bl in 2:16) {
  for (i in seq_len(length(all_net_result_complex_[[bl]]$crhs))) {
    if(length(all_net_result_complex_[[bl]]$crhs[[i]])>1){
      nb_crhs = nb_crhs + 1
    }
  }
}
nb_crhs

########### Sauvegarde des données ###########
save(all_net_result_complex, file = "rdata/all_net_result_complex.rda")

########### Sauvegarde des résulats sur les CRHs complexes fusionnés ###############
save(all_net_result_complex_, file = "rdata/all_net_result_complex_.rda")

##### Les CRHs avant degenerescance #########

crhs_inspection_ = matrix(
  0,
  nrow = 7440,
  ncol = 3
)

l = 1

for (bl in 2:16) {
  for (i in seq_len(length(all_net_result_complex[[bl]]$crhs))) {
    if(length(all_net_result_complex[[bl]]$crhs[[i]])>1){
      crhs_inspection_[l, 1] =  dim(all_net_result_complex[[bl]]$crhs[[i]]$mat_incidence)[1]
      crhs_inspection_[l, 2] =  dim(all_net_result_complex[[bl]]$crhs[[i]]$mat_incidence)[2]
      crhs_inspection_[l, 3] =  paste0("Block ", bl, "- CRH ", i)
      l = l +  1
    }
  }
}

# Premier gros CRH
crhs_inspection_[which.max(as.numeric(crhs_inspection_[, 1])), ]

# Deuxième gros CRH
crhs_inspection_[which.max(as.numeric(crhs_inspection_[, 2])), ]


##### Les CRHs apres degenerescance #########
crhs_inspection = matrix(
  0,
  nrow = 1256,
  ncol = 3
)

l = 1

for (bl in 2:16) {
  for (i in seq_len(length(all_net_result_complex_[[bl]]$crhs))) {
    if(length(all_net_result_complex_[[bl]]$crhs[[i]])>1){
      crhs_inspection[l, 1] =  dim(all_net_result_complex_[[bl]]$crhs[[i]]$mat_incidence)[1]
      crhs_inspection[l, 2] =  dim(all_net_result_complex_[[bl]]$crhs[[i]]$mat_incidence)[2]
      crhs_inspection[l, 3] =  paste0("Block ", bl, "- CRH ", i)
      l = l +  1
    }
  }
}

# Premier gros CRH
crhs_inspection[which.max(as.numeric(crhs_inspection[, 1])), ]

# Deuxième gros CRH
crhs_inspection[which.max(as.numeric(crhs_inspection[, 2])), ]

