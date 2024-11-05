## Cette fonction permet de faire la comparaison entre les CRHs et donne les statistiques

compute_comparaison <- function(all_net_result, clu_chrs_result, make_degeneration = FALSE){
  # print(format(Sys.time(), "%H:%M:%S"))
  
  ################ Création de la matrice servqnt à accuaillir les données de la comparaison ########
  
  matLines = 0
  row_names = ""
  blocs = 2:length(all_net_result)
  for (l in blocs) {
    # print(length(all_net_result[[l]]$crhs))
    for (c in seq_len(length(all_net_result[[l]]$crhs))) {
      if(sum(all_net_result[[l]]$crhs[[c]]$mat_incidence) != -1){
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
  for (clus in seq_along(clu_chrs_result)) {
    
    # On parcours les CRHs de chaque cluster dans cette boucle
    for (crh in seq_along(clu_chrs_result[[clus]])) {
      col = col + 1
      mat = clu_chrs_result[[clus]][[crh]]$mat_incidence
      # On reconstruit le réseau afin de recuperer après les arretes et noeuds à des fin de comparaison
      # net_bip_clus <- graph_from_incidence_matrix(mat)
      
      ln = 0
      # On boucle sur les blocks allant de 2 à 16
      for (bl in blocs) {
        # On recupere le CRHs du block courant
        block <- all_net_result[[bl]]
        crhs_count <- length(block$crhs)
        # On boucle sur les CRHs recuperés
        for (crh_ in seq_len(crhs_count)) {
          crh <- block$crhs[[crh_]]
          if(sum(crh$mat_incidence) != -1){
            ln = ln + 1
            mat_degeneration <- crh$mat_incidence
            # Construction de la matrice de degenerescence
            if(make_degeneration){
              mat_degeneration <- degenerationMatrix(crh$mat_incidence)
            }
            
            rowmat = intersect(rownames(mat), rownames(mat_degeneration))
            colmat = intersect(colnames(mat), colnames(mat_degeneration))
            
            # if(length(rowmat)== 0 || length(colmat)== 0){
            #   crhs_comparation_res$sensibility_mat[ln, col] = NA
            #   crhs_comparation_res$specificity_mat[ln, col] = NA
            # }else{
            #   crhs_comparation_res$sensibility_mat[ln, col] = sum(mat_degeneration[rowmat, colmat][mat[rowmat, colmat]==mat_degeneration[rowmat, colmat]]==1)/sum(mat_degeneration==1)
            #   
            #   crhs_comparation_res$specificity_mat[ln, col] = sum(mat_degeneration[rowmat, colmat][mat[rowmat, colmat]==mat_degeneration[rowmat, colmat]]==0)/sum(mat_degeneration==0)
            # }
            
            crhs_comparation_res$sensibility_mat[ln, col] = sum(mat_degeneration[rowmat, colmat][mat[rowmat, colmat]==mat_degeneration[rowmat, colmat]]==1)/sum(mat_degeneration==1)
            
            crhs_comparation_res$specificity_mat[ln, col] = sum(mat_degeneration[rowmat, colmat][mat[rowmat, colmat]==mat_degeneration[rowmat, colmat]]==0)/sum(mat_degeneration==0)
          }
        }
      }
    }
    
  }
  
  # print(format(Sys.time(), "%H:%M:%S"))
  
  return(crhs_comparation_res)
  
}








