
# Ce fichier permet de passer des données de Juicer (les sorties de Juicer) aux matrices d'incidence puis stocke les résultats dans le fichier Rda "merge_loops_clus"

merge_loops_clus = list()

######## Création d'une matrice donnant les régions de chaque bac
nb_clusters = 25
mat_dim = 562
start = 109000000
bins = 2667
bins_per_bacs = 4
queryIRanges = matrix(
  0,
  ncol = 2,
  nrow = mat_dim
)

for (bac in seq_len(mat_dim)) {
  
  bac_start = start + (bac-1)*bins*bins_per_bacs
  bac_end = start + bac*bins*bins_per_bacs
  
  queryIRanges[bac, 1] = bac_start
  queryIRanges[bac, 2] = bac_end
}

rangeBins_qry = IRanges(start = queryIRanges[, 1], end = queryIRanges[, 2])
########################## Trouver les overlap et construire les matrices ############################

for (clus in seq_len(nb_clusters)) {
  
  data <- read_table(paste0("rdata/juicerOutputFiles/cluster_", clus,"_hiccups_loops/merged_loops.bedpe"))
  data <- data[data$`#chr1`!="#", ]
  data$x1 <- as.numeric(data$x1)
  data$x2 <- as.numeric(data$x2)
  data$y1 <- as.numeric(data$y1)
  data$y2 <- as.numeric(data$y2)
  
  rangeBins1_sbj = IRanges(start=data$x1, end=data$x2)
  rangeBins2_sbj = IRanges(start=data$y1, end=data$y2)
  
  hts1 = findOverlaps(rangeBins_qry, rangeBins1_sbj)
  hts1 = as.data.frame(hts1)
  hts2 = findOverlaps(rangeBins_qry, rangeBins2_sbj)
  hts2 = as.data.frame(hts2)
  
  juicer_hts_regions = unique(hts1$subjectHits)
  
  mat_dim = 562
  incidence_mat = matrix(
    0,
    nrow = mat_dim,
    ncol = mat_dim,
    dimnames = list(
      paste0("BAC", sprintf("%03d", seq_len(mat_dim))), 
      paste0("BAC", sprintf("%03d", seq_len(mat_dim)))
    )
  )
  
  for (reg in juicer_hts_regions) {
    
    reg_x = hts1[hts1$subjectHits==reg, ]$queryHits
    reg_y = hts2[hts2$subjectHits==reg, ]$queryHits
    if(length(reg_x) != 0 & length(reg_y) !=0){
      incidence_mat[reg_x, reg_y] = 1
      incidence_mat[reg_y, reg_x] = 1
    }
    
  }
  
  merge_loops_clus[[length(merge_loops_clus)+1]] = incidence_mat
  names(merge_loops_clus) <- paste0("cluster_", 1:length(merge_loops_clus))
  
}

save(merge_loops_clus, file = paste0("rdata/all_rda_data/merge_loops_clus_", nb_clusters, ".rda"))
