# Chargelent des données

load("rdata/clu_chrs_result.rda")

# Inspection des CRHS
ncrhs = 17

crhs_inspections = matrix(
  0,
  nrow = ncrhs,
  ncol = 4
)

ln = 1
for (crh in seq_len(length(clu_chrs_result))) {
  
  for (crh_ in seq_len(length(clu_chrs_result[[crh]]$crhs))) {
    
    crhs_inspections[ln, 1] = nrow(clu_chrs_result[[crh]]$crhs[[crh_]]$mat_incidence)
    crhs_inspections[ln, 2] = ncol(clu_chrs_result[[crh]]$crhs[[crh_]]$mat_incidence)
    crhs_inspections[ln, 3] = sum(clu_chrs_result[[crh]]$crhs[[crh_]]$mat_incidence)
    if((crhs_inspections[ln, 1] == 1 & crhs_inspections[ln, 2] == 1) | (crhs_inspections[ln, 1] == 2 & crhs_inspections[ln, 2] == 1) | (crhs_inspections[ln, 1] == 1 & crhs_inspections[ln, 2] == 2)){
      crhs_inspections[ln, 4] = 1
    }
    ln = ln + 1
  }
  
}

crhs_inspections

# Combien de CRHs moins complexes ?
sum(crhs_inspections[, 4]==1)



