

# Chargement des données
source("./scripts/load_save_data.R")

# Recupération du nombre de structure
nb_replicas = length(unique(all_paired_structure$paire))

nb_replicas = 1
# La boucle ira de 1 à 500
#nb_replicas = 50
for (r in 1:nb_replicas) {
  
  # Construction de la matrice d'incidence avec la première structure de la cellule
  str = 1
  data <- all_paired_structure%>%
    filter(paire == str_c(r, sprintf("%03d", str)))%>%
    select(-c(ends_with("_c"), "paire"))%>%
    mutate(
      ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n()))
    )%>%
    filter(X4 != 1)
  
  three_dim_data <- data[, 1:3]
  mat_row_name = data[, 5]
  compute_dist <- as.matrix(
    dist(
      x = three_dim_data, 
      method = "euclidean", 
      diag = TRUE
    )
  )
  
  compute_dist_bin1 <- sapply(
    compute_dist, 
    FUN = function(item) ifelse(item<3, 1, 0)
  )
  
  compute_dist_bin1 <- matrix(
    compute_dist_bin1, 
    nrow = nrow(compute_dist), 
    ncol = ncol(compute_dist), 
    byrow = TRUE,
    dimnames = list(mat_row_name, mat_row_name)
  )
  
  diag(compute_dist_bin1) <- 0
  
  # compute_dist_bin1 est la première matrice d'incidence. Il faut mettre à présent à 0 toutes les interaction entre les billes de différentes couleures de la matrice. Il nous suffit donc de verifier si les billes en ligne et en colonne sont du meme bloc
  for (l in mat_row_name) {
    for (c in mat_row_name) {
      if (str_sub(l, 1, 3) != str_sub(c, 1, 3)){
        compute_dist_bin1[l, c] = 0
      }
    }
  }
}

# Nombre de contact trouvé dans la matrice d'incidence
sum(compute_dist_bin1)
# La matrice est bien symétrique
isSymmetric(compute_dist_bin1)

# Sauvegarde le fichier crée
save(compute_dist_bin1, file = "rdata/compute_dist_bin1_cell1_.rda")

load("rdata/compute_dist_bin1_cell1_.rda")

stat <- apply(compute_dist_bin1, 1, sum)

# La proportion du nombre de contact recupéré
2706/(sum(compute_dist_bin1)*2)


# Le nombre de contacts dans les 50 matrices simulées
matrix_sum_scenario0 = rep(0, 50)
for (i in 1:50) {
  matrix_sum_scenario0[i] = sum(
    read.table(paste0("Sc_Cell_In-silico_Hi-C/scenario0/hic_mat_", sprintf("%03d", i), ".txt"), quote="\"", comment.char="", stringsAsFactors = FALSE)
  )
}

summary(matrix_sum_scenario0)

summary(matrix_sum_scenario0/(16738*2))


######################################## compute incidence matrix per bac

incidence_matrice_fn <- function(structure, beads_to_bins){
  beads_number = nrow(structure)
  bins_number = beads_number%/%beads_to_bins
  
  matrix_incidence <- matrix(
    0, 
    nrow = bins_number, 
    ncol = bins_number, 
    byrow = TRUE
  )
  structure$ID <- paste0("B", sprintf("%02d", structure$V4), sprintf("%04d", 1:beads_number))
  
  if(beads_number%%beads_to_bins==0){
    
    for (bins_i in seq(bins_number)) {
      start_i = bins_i*beads_to_bins - beads_to_bins + 1
      end_i = bins_i*beads_to_bins
      
      for (bins_j in seq(bins_number)) {
        start_j = bins_j*beads_to_bins - beads_to_bins + 1
        end_j = bins_j*beads_to_bins
        
        if(bins_i == bins_j){
          matrix_incidence[bins_i, bins_j] = 0
        }else{
          data_in_bins = bind_rows(
            structure[start_i:end_i, ], 
            structure[start_j:end_j, ]
          )
          
          three_dim_data <- data_in_bins[, 1:3]
          mat_row_name = data_in_bins$ID
          
          compute_dist <- as.matrix(
            dist(
              x = three_dim_data, 
              method = "euclidean", 
              diag = TRUE
            )
          )
          
          compute_dist_bin <- sapply(
            compute_dist, 
            FUN = function(item) ifelse(item<3, 1, 0)
          )
          
          compute_dist_bin <- matrix(
            compute_dist_bin, 
            nrow = nrow(compute_dist), 
            ncol = ncol(compute_dist), 
            byrow = TRUE,
            dimnames = list(mat_row_name, mat_row_name)
          )
          
          diag(compute_dist_bin) <- 0
          # On retire les contacts entre les billes d'un meme bac
          compute_dist_bin[1:beads_to_bins, 1:beads_to_bins] <- 0
          bac2_i = beads_to_bins + 1
          bac2_j = beads_to_bins*2
          compute_dist_bin[bac2_i:bac2_j, bac2_i:bac2_j] <- 0
          
          for (l in mat_row_name) {
            for (c in mat_row_name) {
              # On retire les contacts entre les billes du premiers block puisqu'elles sont inertes
              if(str_sub(l, 1, 3) == "B01"){
                compute_dist_bin[l, c] <- 0
              }
              if(str_sub(c, 1, 3) == "B01"){
                compute_dist_bin[l, c] <- 0
              }
              # On retire les contacts des billes qui ne sont pas dans le meme block
              if (str_sub(l, 1, 3) != str_sub(c, 1, 3)){
                compute_dist_bin[l, c] <- 0
              }
            }
          }
          matrix_incidence[bins_i, bins_j] <- sum(compute_dist_bin)
          matrix_incidence[bins_j, bins_i] <- sum(compute_dist_bin)
        }
      }
      
    }
    return(matrix_incidence)
  }else{
    print("le nombre de bacs n'est pas divisible avec le nombre de billes !")
  }
}






incidence_matrice_fn2 <- function(structure, beads_to_bins){
  beads_number = nrow(structure)
  bins_number = beads_number%/%beads_to_bins
  
  matrix_incidence <- matrix(
    0, 
    nrow = bins_number, 
    ncol = bins_number, 
    byrow = TRUE
  )
  structure$ID <- paste0("B", sprintf("%02d", structure$V4), sprintf("%04d", 1:beads_number))
  
  if(beads_number%%beads_to_bins==0){
    
    for (bins_i in seq(bins_number)) {
      start_i = bins_i*beads_to_bins - beads_to_bins + 1
      end_i = bins_i*beads_to_bins
      
      for (bins_j in seq(bins_number)) {
        start_j = bins_j*beads_to_bins - beads_to_bins + 1
        end_j = bins_j*beads_to_bins
        
        if(bins_i != bins_j){
          data_in_bins = bind_rows(
            structure[start_i:end_i, ], 
            structure[start_j:end_j, ]
          )
          mat_row_name = data_in_bins$ID
          
          for (beads_i in seq_len(beads_to_bins)) {
            for (beads_j in seq_len(beads_to_bins)) {
              # On ne retient pas les contacts entre les billes du premiers block puisqu'elles sont inertes. Donc on s'assure qu'elles ne viennent pas du block 1
              # On ne retient pas les contacts des billes qui ne sont pas dans le meme block. Donc on s'assure que les billes des deux bacs soient dans le meme block
              if(str_sub(mat_row_name[beads_i], 1, 3) != "B01" & str_sub(mat_row_name[beads_j+4], 1, 3) != "B01" & str_sub(mat_row_name[beads_i], 1, 3) != str_sub(mat_row_name[beads_j+4], 1, 3) & dist(data_in_bins[c(beads_i, beads_j+4), 1:3], method = "euclidean") < 3){
                matrix_incidence[bins_i, bins_j] = matrix_incidence[bins_i, bins_j] + 1
                matrix_incidence[bins_j, bins_i] = matrix_incidence[bins_j, bins_i] + 1
                
              }
            }
          }
        }
      }
      
    }
    return(matrix_incidence)
  }else{
    print("le nombre de bacs n'est pas divisible avec le nombre de billes !")
  }
}


mat <- incidence_matrice_fn(structure, 4)

mat <- incidence_matrice_fn1(data, 4)

for (i in 1:10) {
  data <- read.table(paste0("hic_mat_", sprintf("%03d", i), ".txt"), quote="\"", comment.char="", stringsAsFactors = FALSE)
  print(max(data))
  print(sum(data))
  print(dim(data))
}




























