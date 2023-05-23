# Disposition des données dans une liste
ncol = 562
ncells = 250
data_dim = ncol*(ncol-1)/2

cellUperDiagData <- matrix(
  0,
  nrow = data_dim,
  ncol = ncells,
  byrow = TRUE
)

for (i in seq_len(ncells)) {
  mat <- as.matrix(
    read.table(
      paste0("hic_mat_cell_", sprintf("%03d", i), ".txt"),
      quote="\"",
      comment.char="",
      stringsAsFactors = FALSE
    )
  )
  
  cellUperDiagData[, i] <- mat[upper.tri(mat)]
}
rm("mat")

save(cellUperDiagData_with_rep, file = "cellUperDiagData_with_rep.rda")

