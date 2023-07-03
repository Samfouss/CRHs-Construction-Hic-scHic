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
      paste0("rdata/single_cell_hic_data/hic_ideal_data_with_rep/hic_mat_cell_", sprintf("%03d", i), ".txt"),
      quote="\"",
      comment.char="",
      stringsAsFactors = FALSE
    )
  )
  
  cellUperDiagData[, i] <- mat[upper.tri(mat)]
}
rm("mat")

save(cellUperDiagData_with_rep, file = "rdata/all_rda_data/cellUperDiagData_with_rep.rda")


################################ Données non ideal à transformer pour le clustering
#library(bigmemory)
ncol = 562
ncells = 250
data_dim = ncol*(ncol-1)/2
nrep = 40

cellUperDiagData_with_repT1 <- matrix(
  0,
  nrow = (ncells*nrep)/2,
  ncol = data_dim
)
l = 1

cells = seq_len(ncells)
repts = seq_len(nrep)

for (i in cells[1:length(cells)]) {
  for( j in repts[1:(length(repts)/2)]){
    mat <- as.matrix(
      read.table(
        #paste0("rdata/single_cell_hic_data/hic_data_with_rep/hic_mat_cell_", sprintf("%03d", i),"_", sprintf("%02d", j),".txt"),
        paste0("../BackUp14022020/rdata/single_cell_hic_data_with_rep/", "hic_mat_cell_", sprintf("%03d", i),"_", sprintf("%02d", j),".txt"),
        quote="\"",
        comment.char="",
        stringsAsFactors = FALSE
      )
    )
    
    cellUperDiagData_with_repT1[l, ] <- mat[upper.tri(mat)]
    rm("mat")
    l = l + 1
  }
}

cellUperDiagData_with_repT2 <- matrix(
  0,
  nrow = (ncells*nrep)/2,
  ncol = data_dim
)

for (i in cells[1:length(cells)]) {
  for( j in repts[(length(repts)/2 + 1):length(repts)]){
    mat <- as.matrix(
      read.table(
        #paste0("rdata/single_cell_hic_data/hic_data_with_rep/hic_mat_cell_", sprintf("%03d", i),"_", sprintf("%02d", j),".txt"),
        paste0("../BackUp14022020/rdata/single_cell_hic_data_with_rep/", "hic_mat_cell_", sprintf("%03d", i),"_", sprintf("%02d", j),".txt"),
        quote="\"",
        comment.char="",
        stringsAsFactors = FALSE
      )
    )
    
    cellUperDiagData_with_repT2[l, ] <- mat[upper.tri(mat)]
    rm("mat")
    l = l + 1
  }
}

save(cellUperDiagData_with_repT, file = "rdata/all_rda_data/cellUperDiagData_with_repT.rda")
