library(tidyverse)
library(data.table)
####### Chargement des données et sauvegarde pour une utilisation ultérieure

# On charge ici d'abord les données sur les structures non pairées (singletons)
all_single_structure <- data.frame()
for (j in 1:500) {
  data_url = paste0("https://raw.githubusercontent.com/fmusella/In-silico_Hi-C_GAM_SPRITE/main/data/structures/structure_", j, ".txt")
  
  all_single_structure <- bind_rows(
    all_single_structure,
    read_table(
      data_url, 
      col_names = FALSE
    )%>%
      mutate(
        replica = sprintf("%03d", j),
        X1_c = X1 - mean(X1),
        X2_c = X2 - mean(X2),
        X3_c = X3 - mean(X3)
      )
  )
}

# On charge ici d'abord les données sur les structures pairées
all_paired_structure <- data.frame()
for (j in 1:250) {
  data_url = paste0("https://raw.githubusercontent.com/fmusella/In-silico_Hi-C_GAM_SPRITE/main/data/pairs/pair_", j, ".txt")
  data_name = paste("structure_", j, sep="")
  
  all_paired_structure <- bind_rows(
    all_paired_structure,
    read_table(data_url, col_names = FALSE)%>%
      select(c("X1", "X2", "X3", "X4"))%>%
      mutate(
        paire = paste0(1, sprintf("%03d", j)),
        X1_c = X1 - mean(X1),
        X2_c = X2 - mean(X2),
        X3_c = X3 - mean(X3)
      ),
    read_table(data_url, col_names = FALSE)%>%
      select(c("X5", "X6", "X7", "X8"))%>%
      rename(X1 = X5, X2 = X6, X3 = X7, X4 = X8)%>%
      mutate(
        paire = paste0(2, sprintf("%03d", j)),
        X1_c = X1 - mean(X1),
        X2_c = X2 - mean(X2),
        X3_c = X3 - mean(X3)
      )
  )
}

for (j in 1:250) {
  data_url = paste0("https://raw.githubusercontent.com/fmusella/In-silico_Hi-C_GAM_SPRITE/main/data/pairs/pair_", j, ".txt")
  data_name = paste("structure_", j, sep="")
  
  folder = str_c("cell_", sprintf("%03d", j))
  paire_1 = read_table(data_url, col_names = FALSE)%>%
    select(c("X1", "X2", "X3", "X4"))
  
  paire_2 = read_table(data_url, col_names = FALSE)%>%
    select(c("X5", "X6", "X7", "X8"))%>%
    rename(X1 = X5, X2 = X6, X3 = X7, X4 = X8)
  
  if(file.exists(str_c("rdata/", folder))){
    
    write.table(
      x = paire_1, 
      file = str_c("rdata/", folder, "/paire_1.txt"), 
      fileEncoding = "UTF-8",
      col.names=FALSE,
      row.names=FALSE
    )
    write.table(
      x = paire_2, 
      file = str_c("rdata/", folder, "/paire_2.txt"), 
      fileEncoding = "UTF-8",
      col.names=FALSE,
      row.names=FALSE
    )
    
  }else{
    
    dir.create(str_c("rdata/", folder))
    write.table(
      x = paire_1, 
      file = str_c("rdata/", folder, "/paire_1.txt"), 
      fileEncoding = "UTF-8",
      col.names=FALSE,
      row.names=FALSE
    )
    write.table(
      x = paire_2, 
      file = str_c("rdata/", folder, "/paire_2.txt"), 
      fileEncoding = "UTF-8",
      col.names=FALSE,
      row.names=FALSE
    )
  }
  
}

# On sauve les bases de données afin de ne plus avoir à les retelecharger
save(all_paired_structure, file = "rdata/all_paired_structure.rda")
save(all_single_structure, file = "rdata/all_single_structure.rda")

####### Pour trouver les replicats qui matchent, on utilise les étapes suivantes :

# 1- Il faut mettre toutes les données dans un seul data frame en prenant soins d'identifier les replicas
# 2 - Concatener toutes les données (les coordonnées puis la couleur d'appartenance) et ensuite les fusionner en utilisant comme clée la variable contenant les données concatenées
# 3 - Pour les données pairées, on procède comme précédemment.

## Ici on charge les données qu'on a pris soins de charger et souvées au préalable à partir du fichier load_save_data.R
## all_paired_structure.rda contient les données pairées
## all_single_structure.rda contient les données sur les structures
load("rdata/all_paired_structure.rda")
load("rdata/all_single_structure.rda")

### En explorant les données, nous avons constaté que les données issues des structures étaient centrées. Pour trouver lea matches, il fallait ramener les données à la meme mésure. Les données issues des paires ont été donc de ce fait centrées.

# ID conteint les données concatenées, replica le numero du replica des structures non pairées, paire donne l'information sur le numéro des données pairées puis le numéro de la paire.
match <- all_single_structure %>% 
  mutate(
    ID = paste0(round(X1_c, 5), round(X2_c, 5), round(X3_c, 5), X4)
  )%>%
  select(
    c("replica", "ID")
  )%>%
  left_join(
    all_paired_structure%>%
      mutate(
        ID = paste0(round(X1_c, 5), round(X2_c, 5), round(X3_c, 5), X4)
      )%>%
      select(
        c("ID", "paire")
      ),
    by = "ID"
  )

# On compte ici le nombre de points des structures pairées qui n'ont pas pu être appariés : 223
sum(is.na(match$paire))

# L'idée ici est de voir les structures qui s'apparient le plus. perc calcul le pourcentage d'appariement des structures
repl_match <- match%>%
  group_by(replica, paire)%>%
  summarise(perc = (n()/2250)*100)%>%
  arrange(replica)%>%
  filter(!is.na(paire))

View(repl_match)
# On récupère les appariement unique et majoritaire : on a en tout 470 structures qui s'apparient de façon univoque
most_match <- repl_match%>%
  filter(
    !(duplicated(replica) | duplicated(replica, fromLast = TRUE))
  )%>%
  filter(
    perc > 99
  )

# On recupère ici les structures qui ont un taux d'appariement supperieur à 99
others <- repl_match%>%
  filter(
    duplicated(replica) | duplicated(replica, fromLast = TRUE)
  )

others$dup1 <- duplicated(others$paire)
others$dup2 <- duplicated(others$replica)

view(others%>%arrange(replica, paire))
others1 <- others%>%
  filter(dup1 == FALSE)%>%
  filter(dup2 == FALSE)

others2 <- others%>%
  filter(!(replica %in% others1$replica) & !(paire %in% others1$paire))


match <- most_match%>%
  bind_rows(
    others1,
    others2
  )%>%
  select(-c("dup1", "dup2"))%>%
  mutate(
    paired_number = str_sub(paire, 1, 1),
    replica_paired = str_sub(paire, 2, 4)
  )

save(match, file = "rdata/structures_match.rda")

view(match%>%arrange(replica_paired))

load("rdata/structures_match.rda")

structure_data = list()
for (j in 1:250) {
  # Selectionner les lignes des structures pairées dans la base de données match
  slice_str = match%>%
    filter(
      replica_paired == sprintf("%03d", j)
    )
  
  structure_data[[length(structure_data)+1]] <- list(
    "structure_1" = as_tibble(all_single_structure%>%filter(replica == pull(slice_str[1, 1]))),
    "structure_2" = as_tibble(all_single_structure%>%filter(replica == pull(slice_str[2, 1])))
  )
  
}

names(structure_data) <- str_c("paired", "_", 1:250)

# Sauve les données dans l'objet suivant
save(structure_data, file = "rdata/structure_data.rda")

