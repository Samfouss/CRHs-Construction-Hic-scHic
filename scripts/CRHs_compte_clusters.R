# Fonction pour calculer le nombre de clusters dans lequel un CRH se retrouve
CRHs_compte_clusters = function(cel,clusters)
{
  length(unique(clusters[cel=="in"]))
}
# Test
CRHs_compte_clusters(cells_data$crhs001,hcpc_cluster$data.clust$clust)

nclusters_CRHs = apply(cells_data,2,CRHs_compte_clusters,clusters=hcpc_cluster$data.clust$clust)
table(nclusters_CRHs)

# Nombre de CRHs différents par cluster
ncl = nlevels(hcpc_cluster$data.clust$clust)
nCRHs_par_cluster = numeric(ncl)
for (i in 1:ncl)
{
  iCRH = apply(cells_data[hcpc_cluster$data.clust$clust==i,],2,function(vec) any(vec=="in"))
  nCRHs_par_cluster[i] = sum(iCRH)
}