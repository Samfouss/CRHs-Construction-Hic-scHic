inciddon = inciddon + rbinom(prod(dim(incidA)),1,pn)

load("rdata/all_rda_data/clu_chrs_result_3Mb.rda")
load("rdata/all_rda_data/all_net_result_3Mb.rda")

mat = clu_chrs_result$cluster_92$crh2$mat_incidence

# all_net_result_ = all_net_result
all_net_result = all_net_result_3Mb
# all_net_result__ = all_net_result
# all_net_result = all_net_result_

sen = c()
spec = c()
desc = c()

for (bl in 2:16) {
  
  for (i in seq_len(length(all_net_result[[bl]]$crhs))) {
    mat_degeneration = all_net_result[[bl]]$crhs[[i]]$mat_incidence
    mat_degeneration = degenerationMatrix(mat_degeneration)
    
    if(sum(mat_degeneration) != -1 & nrow(mat_degeneration) > 2 & ncol(mat_degeneration) > 2){
      rowmat = union(rownames(mat), rownames(mat_degeneration))
      colmat = union(colnames(mat), colnames(mat_degeneration))
      
      mat_degeneration_redim <- matrix(
        -1,
        ncol=length(colmat), 
        nrow=length(rowmat), 
        dimnames=list(rowmat, colmat)
      )
      mat_redim <- mat_degeneration_redim
      
      indxA <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat), colnames(mat), FUN=paste)
      indxB <- outer(rowmat, colmat, FUN=paste) %in% outer(rownames(mat_degeneration), colnames(mat_degeneration), FUN=paste)
      mat_redim[indxA] <- mat
      mat_degeneration_redim[indxB] <- mat_degeneration
      
      a = sum(mat_degeneration_redim[mat_redim==mat_degeneration_redim]==1)
      b = sum(mat_degeneration_redim==1)
      sen = c(sen, a/b)
      
      
      a_ = sum(mat_degeneration_redim[mat_redim==mat_degeneration_redim]==0)
      b_ = sum(mat_degeneration_redim==0) 
      spec = c(spec, a_/b_)
      if(is.na(a/b)){
        print(paste0("sens : ", a/b, "-", bl, "-", i))
      }
      if(is.na(a_/b_)){
        print(paste0("spec : ", a_/b_, "-", bl, "-", i))
      }
      
      # desc = c(desc, paste0(bl, "-", i))
    }
  }
  
}

length(sen)
length(spec)

sen[is.na(sen)] = 0
spec[is.na(spec)] = 0
# 358 avec 0.2
# Avec 0.5 : 696
length(which(sen>0.5))
# 203 avec 0.2
# Avec 0.5 : 260
length(which(spec>0.5))

which(sen>0.5)
which(spec>0.5)

num = 1198
bl = as.numeric(str_split_1(desc[num], "-")[1])
crh = as.numeric(str_split_1(desc[num], "-")[2])

mat1 = clu_chrs_result$cluster_92$crh2$mat_incidence
mat2 = degenerationMatrix(all_net_result[[14]]$crhs[[58]]$mat_incidence)

View(mat1)
View(mat2)
# On prend les lignes et colonne de la matrice mat2
c = row.names(mat2)[row.names(mat2) %in% row.names(mat1)]
l = colnames(mat2)[colnames(mat2) %in% colnames(mat1)]
View(mat1[c, l])

#############################
load("rdata/all_rda_data/crhs_comparaison_res.rda")
load("rdata/all_rda_data/all_net_result_complex_3Mb_.rda")
# Inspection des résultats
specif = crhs_comparaison_res$specificity_mat

nb_crhs_in_clus = dim(specif)[2]
nb_crh_in_struc = dim(specif)[1]
line_name = row.names(specif)
specif_NA_mat_res = matrix(0, nrow = length(line_name), ncol = 3)

l1 = 1
for (l in seq_len(nb_crh_in_struc)) {
  line_data = specif[l, ]
  name = line_name[l]
  if(all(is.na(line_data))){
    bl = as.numeric(str_split_1(name, "_")[2])
    crh = as.numeric(str_split_1(name, "_")[4])
    mat = all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence
    if(all(mat==1)){specif_NA_mat_res[l1, 1] = 1}
    specif_NA_mat_res[l1, 2] = nrow(mat)
    specif_NA_mat_res[l1, 3] = ncol(mat)
    l1 = l1 + 1
  }
}

specif_NA_mat_res = specif_NA_mat_res[specif_NA_mat_res[, 1]!=0, ]

save(specif_NA_mat_res, file = "rdata/all_rda_data/specif_NA_mat_res.rda")





