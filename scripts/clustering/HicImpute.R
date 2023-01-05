library("devtools")
## Warning: le package 'devtools' a été compilé avec la version R 4.1.3
## Le chargement a nécessité le package : usethis
## Warning: le package 'usethis' a été compilé avec la version R 4.1.3
# Install "HiCImpute" package from github.
# install_github("https://github.com/sl-lin/HiCImpute")
library("HiCImpute")

# Disposition des données dans une liste
ncol = 562
ncells = 250
data_dim = ncol*(ncol-1)/2

hicImpute_data <- matrix(
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
  
  hicImpute_data[, i] <- mat[upper.tri(mat)]
}
rm("mat")


scHic = hicImpute_data[, 1:50]
# HicImpute
MCMCImpute_result=MCMCImpute(
  scHiC=scHic,
  bulk=apply(scHic,1,sum),
  expected=NULL,
  startval=c(100,100,10,8,10,0.1,900,0.2,0,replicate(dim(scHic)[2],8)),
  n=562,
  mc.cores = 1,
  cutoff=0.1,
  niter=1,
  burnin=1
)


data("K562_T1_4k")
max(K562_T1_4k)
table(K562_T1_4k)
25759/(54592+45122+25022+13396+7686+4602+3108+1856+970+450+219+121+60+30+7)

max(hicImpute_data)
table(hicImpute_data)
39111990/(39111990+271279+24491+2274+207+8+1)