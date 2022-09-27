# Fonction pour calculer le nombre de clusters dans lequel un CRH se retrouve

CRHs_compte_clusters = function(cel,clusters)
{
  length(unique(clusters[cel=="in"]))
}
# Fonction pour calculer le nombre de clusters dans lequel se retrouve un CRH présent dans au moins tm cellules 
CRHs_compte_clusters_tm = function(cel,clusters,tm)
{
  sum(table(clusters[cel=="in"])>=tm)
}

# Test
<<<<<<< HEAD
# CRHs_compte_clusters(cells_data$crhs001,hcpc_cluster$data.clust$clust)

# Fonction qui donne le nombre de CRHs différents par cluster

CRHs_unique_compte_clusters <- function(cells_data, hcpc_cluster){
  
  ncl = nlevels(hcpc_cluster$data.clust$clust)
  nCRHs_par_cluster = numeric(ncl)
  
  for (i in 1:ncl)
  {
    iCRH = apply(
      cells_data[hcpc_cluster$data.clust$clust==i,],
      2,
      function(vec) any(vec=="in")
    )
    nCRHs_par_cluster[i] = sum(iCRH)
  }
  
  return(
    nCRHs_par_cluster
  )
}


## Cette fonction utilise les précédentes fonctions afin de retourner le nombre de clusters dans lequel un CRH se retrouve, le nombre de CRHs différents par cluster puis le nombre de cellules par cluster
resume_cluster <- function(cells_data_to_use, hcpc_cluster_to_use){
  
  nclusters_CRHs = apply(
    cells_data_to_use,
    2,
    CRHs_compte_clusters,
    clusters=hcpc_cluster_to_use$data.clust$clust
  )
  res1 <- table(nclusters_CRHs)
  
  res2 <- data.frame(table(hcpc_cluster_to_use$data.clust$clust))
  names(res2)[names(res2) == "Var1"] <- "cluster"
  names(res2)[names(res2) == "Freq"] <- "number_of_cells"
  res2$number_of_unique_crhs <- CRHs_unique_compte_clusters(
    cells_data_to_use, 
    hcpc_cluster_to_use
  )
  
  return(
    list(
      "nclusters_CRHs" = res1,
      "CRHs_cells_compte_clusters" = res2
    )
  )
  
}


test = resume_cluster(cells_data, hcpc_cluster)



||||||| e5b9ec8
CRHs_compte_clusters(cells_data$crhs001,hcpc_cluster$data.clust$clust)
=======
CRHs_compte_clusters(cells_data[,1],hcpc_cluster$data.clust$clust)
CRHs_compte_clusters_tm(cells_data[,1],hcpc_cluster$data.clust$clust,1)
>>>>>>> c0271fd11d9fc3f4b8569922fe331fa1f97f6fb0

<<<<<<< HEAD
||||||| e5b9ec8
nclusters_CRHs = apply(cells_data,2,CRHs_compte_clusters,clusters=hcpc_cluster$data.clust$clust)
table(nclusters_CRHs)
=======
nclusters_CRHs = apply(cells_data,2,CRHs_compte_clusters,clusters=hcpc_cluster$data.clust$clust)
table(nclusters_CRHs)
summary(nclusters_CRHs)

nclusters_CRHs = apply(cells_data_,2,CRHs_compte_clusters,clusters=hcpc_cluster_$data.clust$clust)
table(nclusters_CRHs)

nclusters_CRHst5 = apply(cells_data,2,CRHs_compte_clusters_tm,clusters=hcpc_cluster$data.clust$clust,tm=5)
table(nclusters_CRHst5)
summary(nclusters_CRHst5[nclusters_CRHst5>0])
sum(nclusters_CRHst5>0)
>>>>>>> c0271fd11d9fc3f4b8569922fe331fa1f97f6fb0

<<<<<<< HEAD
||||||| e5b9ec8
# Nombre de CRHs différents par cluster
ncl = nlevels(hcpc_cluster$data.clust$clust)
nCRHs_par_cluster = numeric(ncl)
for (i in 1:ncl)
{
  iCRH = apply(cells_data[hcpc_cluster$data.clust$clust==i,],2,function(vec) any(vec=="in"))
  nCRHs_par_cluster[i] = sum(iCRH)
}
=======
# Nombre de CRHs différents par cluster
ncl = nlevels(hcpc_cluster$data.clust$clust)
nCRHs_par_cluster = numeric(ncl)
for (i in 1:ncl)
{
  iCRH = apply(cells_data[hcpc_cluster$data.clust$clust==i,],2,function(vec) any(vec=="in"))
  nCRHs_par_cluster[i] = sum(iCRH)
}
nCRHs_par_cluster
summary(nCRHs_par_cluster)

# Nombre de clusters avec au moins tm répétitions de chaque CRH
tm=5
nCRHst5_par_cluster = numeric(ncl)
for (i in 1:ncl)
{
  iCRH = apply(cells_data[hcpc_cluster$data.clust$clust==i,],2,function(vec) sum(vec=="in")>=tm)
  nCRHst5_par_cluster[i] = sum(iCRH)
}
nCRHst5_par_cluster
summary(nCRHst5_par_cluster)
>>>>>>> c0271fd11d9fc3f4b8569922fe331fa1f97f6fb0
