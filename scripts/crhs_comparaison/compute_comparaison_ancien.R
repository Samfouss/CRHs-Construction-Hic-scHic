compute_comparaison <- function(all_net_result, clu_chrs_result, make_degeneration = TRUE){
  
  ################ Création de la matrice servqnt à accuaillir les données de la comparaison ########
  
  ### Noms des lignes de la matrice : les CRHs puis leur block d'appartenance
  matLines = 0
  row_names = ""
  for (l in 2:length(all_net_result)) {
    # print(length(all_net_result[[l]]$crhs))
    for (c in seq_len(length(all_net_result[[l]]$crhs))) {
      if(length(all_net_result[[l]]$crhs[[c]]) == 2){
        matLines = matLines +  1
        row_names = c(row_names, str_c("block_", l, "_crhs_", c))
      }
    }
  }
  row_names = row_names[row_names != ""]
  
  ### Noms des colonnes de la matrice : Les chrs puis leur cluster d'appartenance
  clus_num = c()
  for (cell in names(clu_chrs_result)) {
    n = as.numeric(unlist(strsplit(cell, "[_]"))[2])
    clus_num = c(clus_num, n)
  }
  matCol = 0
  col_names = ""
  for (clus in seq_len(length(clu_chrs_result))) {
    # print(length(clu_chrs_result[[clus]]$crhs))
    l = length(clu_chrs_result[[clus]])
    matCol = matCol + l
    col_names = c(col_names, str_c("cluster_", clus_num[clus], "_crhs_", seq_len(l)))
  }
  col_names = col_names[col_names != ""]
  
  ### Création de la matrice pour accueillir les résulats
  crhs_comparation_res = list(
    "sensibility_mat" = matrix(
      0,
      nrow = matLines,
      ncol = matCol,
      dimnames = list(row_names, col_names)
    ),
    "specificity_mat" = matrix(
      0,
      nrow = matLines,
      ncol = matCol,
      dimnames = list(row_names, col_names)
    )
  )
  
  col = 0
  
  # On parcours les cluster dans cette boucle
  for (clus in seq_len(length(clu_chrs_result))) {
    
    # On parcours les CRHs de chaque cluster dans cette boucle
    for (crh in seq_len(length(clu_chrs_result[[clus]]))) {
      col = col + 1
      mat = clu_chrs_result[[clus]][[crh]]$mat_incidence
      # On reconstruit le réseau afin de recuperer après les arretes et noeuds à des fin de comparaison
      net_bip_clus <- graph_from_incidence_matrix(mat)
      
      ln = 0
      # On boucle sur les blocks allant de 2 à 16
      for (bl in 2:length(all_net_result)) {
        # On recupere le CRHs du block courant
        crhs_ = length(all_net_result[[bl]]$crhs)
        # On boucle sur les CRHs recuperés
        for (crh_ in seq_len(crhs_)) {
          if(sum(all_net_result[[bl]]$crhs[[crh_]]$mat_incidence)!=-1){
            ln = ln + 1
            mat_degeneration <- all_net_result[[bl]]$crhs[[crh_]]$mat_incidence
            # Construction de la matrice de degenerescence
            if(make_degeneration){
              mat_degeneration <- degenerationMatrix(all_net_result[[bl]]$crhs[[crh_]]$mat_incidence)
            }
            # redimenssion des deux matrices afin de calcluer les statistiques sur les intersections
            rowmat = union(rownames(mat), rownames(mat_degeneration))
            colmat = union(colnames(mat), colnames(mat_degeneration))
            
            mat_degeneration_redim <- matrix(
              -1,
              ncol=length(colmat),
              nrow=length(rowmat),
              dimnames=list(rowmat, colmat)
            )
            mat_redim <- mat_degeneration_redim
            
            indxA <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat), colnames(mat), FUN=paste)
            indxB <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat_degeneration), colnames(mat_degeneration), FUN=paste)
            mat_redim[indxA] <- mat
            mat_degeneration_redim[indxB] <- mat_degeneration
            
            crhs_comparation_res$sensibility_mat[ln, col] = sum(mat_degeneration_redim[mat_redim==mat_degeneration_redim]==1)/sum(mat_degeneration_redim==1)
            
            # mismatchs = length(rowmat[rowmat != rowmat_])*length(colmat[colmat != colmat_])
            crhs_comparation_res$specificity_mat[ln, col] = sum(mat_degeneration_redim[mat_redim==mat_degeneration_redim]==0)/sum(mat_degeneration_redim==0)
            
          }
        }
      }
      
    }
  }
  
  crhs_comparation_res
  
}