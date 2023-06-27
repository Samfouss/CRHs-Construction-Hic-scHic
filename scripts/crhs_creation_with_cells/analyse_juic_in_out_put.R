

########### Analyse les inputs de Juicer ##############

nclus = 5
route = "ignorng_data/pamk_input_5/"

for (clus in seq_len(nclus)) {
  
  data = read.table(paste0(route, "cluster_", clus,"_juicer_input.txt"))
  print(paste0("################ Cluster", clus))
  print(summary(data$V9))
  
}
