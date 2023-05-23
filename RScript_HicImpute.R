
################################################## Run HicImpute ############################
library("devtools")
## Warning: le package 'devtools' a été compilé avec la version R 4.1.3
## Le chargement a nécessité le package : usethis
## Warning: le package 'usethis' a été compilé avec la version R 4.1.3
# Install "HiCImpute" package from github.
# install_github("https://github.com/sl-lin/HiCImpute")
library("HiCImpute")

# Chargement des données
load("MCMCImpute_input.rda")
ncol = 562

# HicImpute
MCMCImpute_result=MCMCImpute(
  scHiC=scHic,
  bulk=apply(scHic,1,sum),
  expected=NULL,
  startval=c(100,100,10,2,10,0.1,900,0.2,0,replicate(dim(scHic)[2],2)),
  n=ncol,
  mc.cores = 5,
  cutoff=0.5,
  niter=60000,
  burnin=10000
)

save(MCMCImpute_result, file = "MCMCImpute_result.rda")

################################################## Graphics ############################


load("MCMCImpute_result.rda")
library(ggplot2)
library(dplyr)

data_superposition <- function(){
  
  nb_cells = sample(1:250, 1, replace=FALSE)
  data <- data.frame(
    data_without_imp = MCMCImpute_result$scHiC[, nb_cells[1]],
    data_with_imp = MCMCImpute_result$Impute_SZ[, nb_cells[1]],
    cell = nb_cells[1]
  )
  
  # for (cell in seq_len(length(nb_cells)-1)) {
  #   data <- bind_rows(
  #     data,
  #     data.frame(
  #       data_without_imp = MCMCImpute_result$scHiC[, nb_cells[cell + 1]],
  #       data_with_imp = MCMCImpute_result$Impute_SZ[, nb_cells[cell + 1]],
  #       cells = paste("cell_", nb_cells[cell + 1])
  #     )
  #   )
  # }
  
  duage_points <- ggplot(
    data = data,
    mapping = aes(
      x = data_without_imp, 
      y = data_with_imp
    )
  ) + 
    geom_point(
      colour = "blue", 
      size = 3
    ) +
    geom_abline(intercept = 0, slope = 1)+
    labs(
      title = paste("adjustment of imputed data to data without imputations for cell ", nb_cells),
      x = "Data without imputations",
      y = "Data with imputations"
    )
  duage_points
  
  # duage_points + 
  #   facet_wrap(facets = ~ cells)
  
}

data_superposition()
