library("devtools")
## Warning: le package 'devtools' a été compilé avec la version R 4.1.3
## Le chargement a nécessité le package : usethis
## Warning: le package 'usethis' a été compilé avec la version R 4.1.3
# Install "HiCImpute" package from github.
# install_github("https://github.com/sl-lin/HiCImpute")
library("HiCImpute")

# Disposition des données dans une liste

# Disposition des données dans une liste
ncol = 562
ncells = 250
data_dim = ncol*(ncol-1)/2

scHic <- matrix(
  0,
  nrow = data_dim,
  ncol = ncells,
  byrow = TRUE
)

for (i in seq_len(ncells)) {
  mat <- as.matrix(
    read.table(
      paste0("rdata/single_cell_hic_data/hic_mat_", sprintf("%03d", i), ".txt"),
      quote="\"",
      comment.char="",
      stringsAsFactors = FALSE
    )
  )
  
  scHic[, i] <- mat[upper.tri(mat)]
}
mat = NULL


# HicImpute
MCMCImpute_result=MCMCImpute(
  scHiC=scHic,
  bulk=apply(scHic,1,sum),
  expected=NULL,
  startval=c(100,100,10,8,10,0.1,900,0.2,0,replicate(dim(scHic)[2],8)),
  n=ncol,
  mc.cores = 5,
  cutoff=0.5,
  niter=10,
  burnin=1
)

save(MCMCImpute_result, file = "rdata/MCMCImpute_result.rda")
