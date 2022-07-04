# Test de fusion_comp
overlap_edge <- edge_identity_overlap(structure_1_net_bip, structure_2_net_bip)

overlap_edge$edge_overlap1$chevauche_1_2_perc
overlap_edge$edge_overlap1$chev_edge_comp1
overlap_edge$edge_overlap1$chev_edge_comp2
paires.crh <- overlap_edge$edge_overlap1$chev_edge_comp

# Exécution de la fusion des CRH sur le block 1
adjcomp <- fusion_comp(
  structure_1_net_bip, 
  structure_2_net_bip,
  paires.crh
)

structure_net_comp1 = structure_1_net_bip[[1]]
structure_net_comp2 = structure_2_net_bip[[1]]
dist_bin1 = structure_1_net_bip[[1]]$dist_bin
dist_bin2 = structure_2_net_bip[[1]]$dist_bin
ens1 = adjcomp$block1$ens1
ens2 = adjcomp$block1$ens2

# Résumé des différences
sum(abs(adjcomp$block1$adjcomp[ens1,ens1]-dist_bin1[ens1,ens1]))
sum(abs(adjcomp$block1$adjcomp[ens2,ens2]-dist_bin2[ens2,ens2]))

# Fusion des réplicats 1 et 2
cc.list = list()
# Boucle sur les composantes appariées
for (i in 1:ncol(paires.crh))
{
  cc.list[[i]] = fusion_comp(
    structure_1_net_bip,
    structure_2_net_bip,
    paires.crh[,i]
  )
}
# On ajoute les composantes non appariées
j = length(cc.list)
# du réplicat 1
for(i in (1:structure_net_comp1$no_bip)[!(1:structure_net_comp1$no_bip)%in%paires.crh[,1] & structure_net_comp1$csize_bip>2])
{
  j = j+1
  cc.list[[j]] = structure_1_net_bip[[1]]$dist_bin[
    names(structure_net_comp1$membership_bip==i),
    names(structure_net_comp1$membership_bip==i)
  ]
}
# du réplicat 2
for(i in (1:structure_net_comp2$no_bip)[!(1:structure_net_comp2$no_bip)%in%paires.crh[,2] & structure_net_comp2$csize_bip>2])
{
  j = j+1
  cc.list[[j]] = structure_2_net_bip$dist_bin[
    names(structure_net_comp2$membership_bip==i),
    names(structure_net_comp2$membership_bip==i)
  ]
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

