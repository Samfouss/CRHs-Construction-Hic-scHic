# Chargement des librairies
library("devtools")
## Warning: le package 'devtools' a été compilé avec la version R 4.1.3
## Le chargement a nécessité le package : usethis
## Warning: le package 'usethis' a été compilé avec la version R 4.1.3
# Install "HiCImpute" package from github.
# install_github("https://github.com/sl-lin/HiCImpute")
library("HiCImpute")

# scHiC_Kmeans sur les données avec imputations
load("rdata/MCMCImpute_result.rda")
hicImpute_data <- MCMCImpute_result$Impute_SZ
ncells = 250
colnames(hicImpute_data) <- str_c("cell", sprintf("%03d", seq(ncells)))

max_cluters_hicImput = 12
cluster_with_imp=scHiC_Kmeans(
  hicImpute_data, 
  centers=max_cluters_hicImput, 
  nstart=1, 
  iter.max=1000, 
  seed=1
)

table(cluster_with_imp$cluster)

cluster_with_imp$cluster[cluster_with_imp$cluster==1]

save(cluster_with_imp, file = "rdata/cluster_with_imp.rda")

