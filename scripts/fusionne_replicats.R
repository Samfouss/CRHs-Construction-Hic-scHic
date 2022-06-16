# Test de fusion_comp
overlap_edg <- edge_identity_overlap(
  structure_1_net, 
  structure_2_net
)

overlap_edg$edge_overlap1$chevauche_1_2_perc
overlap_edg$edge_overlap1$chev_edge_comp1
overlap_edg$edge_overlap1$chev_edge_comp2[, 2:1]

tmp = rbind(
  overlap_edg$edge_overlap1$chev_edge_comp1,
  overlap_edg$edge_overlap1$chev_edge_comp2[,2:1]
)

rbind(
  overlap_edg$chev_crh_comp1,
  overlap_edg$chev_crh_comp2[,2:1]
)

paires.crh = tmp[duplicated(tmp),]
adjcomp = fusion_comp(structure_1_1_net_comp, structure_2_1_net_comp,structure_1_1_net$dist_bin, structure_2_1_net$dist_bin,paires.crh[1,])

structure_net_comp1 = structure_1_1_net_comp
structure_net_comp2 = structure_2_1_net_comp
dist_bin1 = structure_1_1_net$dist_bin
dist_bin2 = structure_2_1_net$dist_bin

# Résumé des différences
sum(abs(adjcomp-dist_bin1[ens1,ens1]))
sum(abs(adjcomp[ens2,ens2]-dist_bin2[ens2,ens2]))
sum(abs(dist_bin1[ens2,ens2]-dist_bin2[ens2,ens2]))

# Fusion des réplicats 1 et 2
cc.list = list()
# Boucle sur les composantes appariées
for (i in 1:nrow(paires.crh))
{
  cc.list[[i]] = fusion_comp(structure_1_1_net_comp, structure_2_1_net_comp,structure_1_1_net$dist_bin, structure_2_1_net$dist_bin,paires.crh[i,])
}
# On ajoute les composantes non appariées
j = length(cc.list)
# du réplicat 1
for(i in (1:structure_1_1_net_comp$no)[!(1:structure_1_1_net_comp$no)%in%paires.crh[,1] & structure_1_1_net_comp$csize>2])
{
  j = j+1
  cc.list[[j]] = structure_1_1_net$dist_bin[structure_1_1_net_comp$membership==i,structure_1_1_net_comp$membership==i]
}
# du réplicat 2
for(i in (1:structure_2_1_net_comp$no)[!(1:structure_2_1_net_comp$no)%in%paires.crh[,2] & structure_2_1_net_comp$csize>2])
{
  j = j+1
  cc.list[[j]] = structure_2_1_net$dist_bin[structure_2_1_net_comp$membership==i,structure_2_1_net_comp$membership==i]
}

# Création d'une nouvelle matrice d'adjacence
nr.vec = c(0,cumsum(sapply(cc.list,nrow)))
nr = nr.vec[length(nr.vec)]
adj_fusion = matrix(0,nr,nr)
for (i in 1:length(cc.list)) adj_fusion[(nr.vec[i]+1):nr.vec[i+1],(nr.vec[i]+1):nr.vec[i+1]] = cc.list[[i]]

structure_fusion_1_net = list("dist_bin" = adj_fusion,"net" = graph_from_adjacency_matrix(
  adj_fusion, mode='undirected'))
par(mar = c(1, 1, 1, 1))
plot(structure_fusion_1_net$net, edge.arrow.size=.4,vertex.label=NA)

