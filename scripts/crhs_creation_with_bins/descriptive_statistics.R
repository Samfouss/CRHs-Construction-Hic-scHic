# Chargement des données
load("rdata/all_net_result_complex_.rda")

# Combien de contacts avons nous en moyenne dans les matrices ?

# Inspection des CRHS
ncrhs = 1173

crhs_inspections_ = matrix(
  0,
  nrow = ncrhs,
  ncol = 4
)

ln = 1
for(bl in 2:length(all_net_result_complex_)) {
  
  for (crh in seq_len(length(all_net_result_complex_[[bl]]$crhs))) {
    if(length(all_net_result_complex_[[bl]]$crhs[[crh]])>1){
      # print(paste(bl, " ", crh, " ", length(all_net_result_complex[[bl]]$crhs[[crh]])))
      crhs_inspections_[ln, 1] = nrow(all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence)
      crhs_inspections_[ln, 2] = ncol(all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence)
      crhs_inspections_[ln, 3] = sum(all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence)
      if(max(all_net_result_complex_[[bl]]$crhs[[crh]]$mat_incidence)>1){
        print(paste(bl, " ", crh))
      }
      if((crhs_inspections_[ln, 1] == 1 & crhs_inspections_[ln, 2] == 1) | (crhs_inspections_[ln, 1] == 2 & crhs_inspections_[ln, 2] == 1) | (crhs_inspections_[ln, 1] == 1 & crhs_inspections_[ln, 2] == 2)){
        crhs_inspections_[ln, 4] = 1
      }
      ln = ln + 1
    }
    
  }
  
}

crhs_inspections_

# Combien de CRHs moins complexes ?
sum(crhs_inspections_[, 4]==1)

mean(crhs_inspections_[, 3])
