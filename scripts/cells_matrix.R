library(tidyverse)
library(ggplot2)

load("rdata/all_net_result.rda")
# On efface tout ce qui est espaces, parenthèses, tirets puis on retient juste le premier nombre de chaque couple dans les résultats

# Construction de la matrice brute
cells_matrix = list()

for (i in seq_len(length(all_net_result))) {
  
  result <- all_net_result[[i]]$resume_fusion%>%
    str_split("-", simplify = TRUE)
  
  result <- gsub(",.*","",result)
  result <- gsub(".* ","",result)
  
  crhs_fus_len = nrow(result)
  
  compute_cells_matrix = matrix(
    0,
    nrow = 250,
    ncol = crhs_fus_len,
    dimnames = list(
      str_c("cell", sprintf("%03d", 1:250)),
      str_c("crhs", sprintf("%03d", 1:crhs_fus_len), "_block", sprintf("%02d", i))
    )
  )
  
  for(r in seq_len(crhs_fus_len)){
    
    # Pour une ligne donnée dans la metrice "result" des résultats, j'extrais les cellules conrespondant aux CRHs fusionnés
    extract_cells = str_sub(result[r, ][result[r, ] != ""], 1, 3)
    
    # Comme on ne compte qu'une seule fois les cellules meme si les deux chromosomes entrent dans la fusion, alors il nous suffit d'effacer les duplicate et de construire une matrice d'incidence à ajouter a la matrice finale "compute_cells_matrix"
    extract_cells = extract_cells[!duplicated(extract_cells)]
    
    compute_cells_matrix[
      str_c("cell", extract_cells),
      r
    ] <- 1
  }
  
  cells_matrix[[length(cells_matrix) + 1]] <- compute_cells_matrix
  
}

names(cells_matrix) <- str_c("block", 1:16)


view(cells_matrix$block1)

################### Distribution ################
nrep_CRHs = apply(cells_matrix$block1,2,sum)
hist(nrep_CRHs)
summary(nrep_CRHs)

ggplot(
  data.frame(nrep_CRHs)%>%
    mutate(nrep_CRHs = as.numeric(nrep_CRHs))%>%
    rownames_to_column( var = "CRHs_names")%>% 
    filter(nrep_CRHs>20)
)+
  geom_bar(mapping = aes(x = CRHs_names, y = nrep_CRHs), stat = "identity")


save(cells_matrix, file = "rdata/cells_matrix.rda")
#save(cells_matrix_, file = "rdata/cells_matrix_.rda")



