#' Simple fonction qui calcul la distance (eucludienne) entre deux points
#'
#' @description
#' #' `dist_func` retourne la distance entre deux points
#'
#' @param point_1 
#' @param point_2 

dist_func <- function(point_1, point_2, dist_method="eucludien"){
  if(dist_method=="eucludien"){
    D = sqrt(apply((point_1 - point_2)^2, MARGIN = 1, FUN = sum))
  }
}

#' Cette fonction calcul la distance entre deux CRHs
#'
#' @description
#' #' `m_dist_func` Cette fonction calcul la distance entre deux CRHs, soit par 
#' rapport aux centres de gravité ou soit point par point
#'
#' @param data_frame1 
#' @param data_frame2
#' 

m_dist_func <- function(data_frame1, data_frame2, type = "diff_mean", output = "matrix", rep1 = 1, rep2 = 2){
  n1 = nrow(data_frame1)
  n2 = nrow(data_frame2)
  
  if(output == "matrix"){
    mat.distance.diff = matrix(NA, nrow = n1, ncol = n2)
    if(type == "diff_mean"){
      row_names1 <- paste0("R", rep1, sprintf("%02d", pull(data_frame1[, 4])))
      row_names2 <- paste0("R", rep2, sprintf("%02d", pull(data_frame2[, 4])))
    }else if(type == "diff_ind"){
      row_names1 <- paste0("R", rep1, "B", sprintf("%02d", data_frame1[, 4]), "CRH", sprintf("%02d", data_frame1[, 6]))
      row_names2 <- paste0("R", rep2, "B", sprintf("%02d", data_frame2[, 4]), "CRH", sprintf("%02d", data_frame2[, 6]))
    }
    
    rownames(mat.distance.diff) <- row_names1
    colnames(mat.distance.diff) <- row_names2
    
  }else if(output == "data.frame"){
    mat.distance.diff = matrix(NA, nrow = n1*n2, ncol = 2)
  }
  
  if(output == "matrix"){
    for (i in seq_len(n1)) {
      for (j in seq_len(n2)) {
        mat.distance.diff[i, j] = round(dist_func(data_frame1[i, 1:3], data_frame2[j, 1:3]), 2)
      }
    }
  }else if(output == "data.frame"){
    #### Calcul des différences de distances moyenne par CRH ####
    if(type == "diff_mean"){
      for (i in 1:n1) {
        for (j in 1:n2) {
          mat.distance.diff[n2*i-n2+j, 1] = paste0("R", rep1, sprintf("%02d", i),"-","R", rep2, sprintf("%02d", j))
          mat.distance.diff[n2*i-n2+j, 2] = round(dist_func(data_frame1[i, 1:3], data_frame2[j, 1:3]), 2)
        }
      }
      #### Calcul des différences de distances par bille ####
    }else if(type == "diff_ind"){
      for (i in 1:n1) {
        for (j in 1:n2) {
          mat.distance.diff[n2*i-n2+j, 1] <- paste0(
            "R", rep1, "B", sprintf("%02d", data_frame1[i, 4]), "CRH", sprintf("%02d", data_frame1[i, 6]),
            "-",
            "R", rep2, "B", sprintf("%02d", data_frame2[i, 4]), "CRH", sprintf("%02d", data_frame2[i, 6])
          )
          mat.distance.diff[n2*i-n2+j, 2] <- round(dist_func(data_frame1[i, 1:3], data_frame2[j, 1:3]), 2)
        }
      }
    }
    mat.distance.diff = data.frame(mat.distance.diff)%>%
      rename(
        "IDs" = "X1",
        "diff_distances" = "X2"
      )
    mat.distance.diff[, 2] <- as.numeric(mat.distance.diff[, 2])
  }
  
  mat.distance.diff
}