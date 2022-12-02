library(tidyverse)

nb_replicas = 250

all_paired_structure = as_tibble(
  read.table(
    paste0("rdata/cell_folder/cell_001/paire_", 1,".txt"), 
    quote="\"", 
    comment.char="", 
    stringsAsFactors = FALSE
  )
)%>%
  mutate(
    X1 = V1,
    X2 = V2,
    X3 = V3,
    X4 = V4,
    ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n())),
    paire = paste0(1, sprintf("%03d", 1))
  )%>%
  select(-starts_with("V"))%>%
  bind_rows(
    as_tibble(
      read.table(
        paste0("rdata/cell_folder/cell_001/paire_", 2,".txt"), 
        quote="\"", 
        comment.char="", 
        stringsAsFactors = FALSE
      )
    )%>%
      mutate(
        X1 = V1,
        X2 = V2,
        X3 = V3,
        X4 = V4,
        ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n())),
        paire = paste0(2, sprintf("%03d", 1))
      )%>%
      select(-starts_with("V"))
  )
  
for (r in 2:nb_replicas) {
  
  # Construction de la matrice d'incidence avec la première structure de la cellule
  for (p in 1:2) {
    all_paired_structure <- all_paired_structure%>%
      bind_rows(
        as_tibble(
          read.table(
            paste0("rdata/cell_folder/cell_", sprintf("%03d", r), "/paire_", p,".txt"), 
            quote="\"", 
            comment.char="", 
            stringsAsFactors = FALSE
          )
        )%>%
          mutate(
            X1 = V1,
            X2 = V2,
            X3 = V3,
            X4 = V4,
            ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n())),
            paire = paste0(p, sprintf("%03d", r))
          )%>%
          select(-starts_with("V"))
    )
  }
  
}

# Sauve les données dans l'objet suivant
save(all_paired_structure, file = "rdata/all_paired_structure.rda")
