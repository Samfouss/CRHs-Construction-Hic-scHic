# Découpage de la régions de 6 Mb en 3 régions de 2 Mb

library(igraph)
load("cellUperDiagData.rda")
# Chargement des promoters
load("/Users/alexandrebureau/Documents/GitHub/CRHs-Construction-Hic-scHic/rdata/all_rda_data/scHic_promoters_ids.rda")

l.entiere = list()
l.sousmat = list()
# Somme de 40 matrices
for (r in 1:250)
{
  cat(r,"\n")
mat = mat_bip = matrix(0,562,562)
mat[upper.tri(mat)] = apply(cellUperDiagData[,(40*(r-1)+1):(40*r)],1,sum)
mat[lower.tri(mat)] = t(mat)[lower.tri(mat)]
dimnames(mat) = list(1:562,1:562)

# Matrice comprenant seulement les connexions entre enhancers et promoteurs
mat_bip = mat[scHic_promoters_ids, -scHic_promoters_ids]

net_mat = graph_from_incidence_matrix(mat_bip)
#plot(net_mat)
net_components_mat <- components(net_mat)
#cat(unique(net_components_mat$csize),"\n")
l.entiere[[r]] = net_components_mat$csize[net_components_mat$csize>1]

# Premiers 2 Mb
prow = max(as.numeric(rownames(mat_bip)[as.numeric(rownames(mat_bip))<=187]))
pcol = max(as.numeric(colnames(mat_bip)[as.numeric(colnames(mat_bip))<=187]))
sousmat1 = mat_bip[1:which(rownames(mat_bip)==prow),1:which(colnames(mat_bip)==pcol)]

# Derniers 2 Mb après bac 375
lrow = min(as.numeric(rownames(mat_bip)[as.numeric(rownames(mat_bip))>=375]))
lcol = min(as.numeric(colnames(mat_bip)[as.numeric(colnames(mat_bip))>=375]))
sousmat2 = mat_bip[(which(rownames(mat_bip)==prow)+1):(which(rownames(mat_bip)==lrow)-1),(which(colnames(mat_bip)==pcol)+1):(which(colnames(mat_bip)==lcol)-1)]
sousmat3 = mat_bip[which(rownames(mat_bip)==lrow):nrow(mat_bip),which(colnames(mat_bip)==lcol):ncol(mat_bip)]

net_mat = graph_from_incidence_matrix(sousmat1)
#plot(net_mat)
net_components_mat <- components(net_mat)
#cat(unique(net_components_mat$csize),"\n")
l1 = unique(net_components_mat$csize)[net_components_mat$csize>1]

net_mat = graph_from_incidence_matrix(sousmat2)
#plot(net_mat)
net_components_mat <- components(net_mat)
#cat(unique(net_components_mat$csize),"\n")
l2 = unique(net_components_mat$csize)[net_components_mat$csize>1]

net_mat = graph_from_incidence_matrix(sousmat3)
#plot(net_mat)
net_components_mat <- components(net_mat)
#cat(unique(net_components_mat$csize),"\n")
l3 = unique(net_components_mat$csize)[net_components_mat$csize>1]
l.sousmat[[r]] = list(l1,l2,l3)
}

sum(sapply(l.entiere,length))
[1] 756
table(sapply(l.entiere,length))
1  2  3  4  5  6  7  8  9 10 13 
41 60 67 50 20  5  1  2  2  1  1 
# Semblable aux résultats de Fousseni
summary(sapply(l.entiere,length))
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
1.000   2.000   3.000   3.024   4.000  13.000 

table(unlist(l.entiere))
summary(unlist(l.entiere))
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
2.0     4.0   131.0   119.2   183.0   402.0 

table(unlist(l.sousmat))
summary(unlist(l.sousmat))
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NAs 
    1.0     7.0    32.0    63.7   137.0   160.0      10 
length(unlist(l.sousmat))
[1] 1336
table(sapply(l.sousmat,function(l) length(l[[1]][l[[1]]>1]) + length(l[[2]][l[[2]]>1]) + length(l[[3]][l[[3]]>1]) ))
3  4  5  6  7  8  9 10 12 14 
17 62 77 62 18  7  3  1  2  1 

