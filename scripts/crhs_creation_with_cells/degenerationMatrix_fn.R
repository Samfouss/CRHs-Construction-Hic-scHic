
degenerationMatrix <- function(matrixToDegenerate){
  
  matrixToDegenerate = matrixToDegenerate
  
  bins_number = 1:2250
  bacs_number = c(0, which(bins_number%%4==0))
  bacs_matrix = matrix(
    "", 
    nrow = length(bins_number), 
    ncol = 2
  )
  
  for (line in seq_len(length(bins_number))) {
    b = which(bacs_number == (line-line%%4))
    if(line%%4==0){
      b = which(bacs_number == (line-line%%4)) - 1
    }
    bacs_matrix[line, 1] = sprintf("%04d", line)
    bacs_matrix[line, 2] = paste0("BAC", sprintf("%03d", b))
  }
  
  row_names = row.names(matrixToDegenerate)
  col_names = colnames(matrixToDegenerate)
  bins_color = unique(str_sub(colnames(matrixToDegenerate), 1, 3))
  row_names_numb = str_sub(row_names, 4, 7)
  col_names_numb = str_sub(col_names, 4, 7)
  
  row_names_ = unique(bacs_matrix[which(bacs_matrix[, 1]%in%row_names_numb), 2])
  col_names_ = unique(bacs_matrix[which(bacs_matrix[, 1]%in%col_names_numb), 2])
  
  mat_reduct <- matrix(
    0, 
    nrow = length(row_names_), 
    ncol = length(col_names_),
    dimnames = list(
      row_names_, 
      col_names_
    )
  )
  
  for (row in row_names_) {
    for (col in col_names_) {
      i = which(row_names %in% paste0(bins_color, bacs_matrix[bacs_matrix[, 2]==row, 1]))
      j = which(col_names %in% paste0(bins_color, bacs_matrix[bacs_matrix[, 2]==col, 1]))
      mat_reduct[row, col] = sum(matrixToDegenerate[i, j])
    }
  }
  
  mat_reduct
  
}
