
get_out_crhs_less_compl <- function(all_net_result){
  
  nb = 0
  for (i in 2:length(all_net_result)) {
    for (k in seq_len(length(all_net_result[[i]]$crhs))) {
      if(length(all_net_result[[i]]$crhs[[k]])>1){
        mat = all_net_result[[i]]$crhs[[k]]$mat_incidence
        # Identification des CRHs moins complex
        if((nrow(mat) == 1 & ncol(mat) == 1) | (nrow(mat) == 2 & ncol(mat) == 1) | (nrow(mat) == 1 & ncol(mat) == 2)){
          all_net_result[[i]]$crhs[[k]] <- -1
          all_net_result[[i]]$resume_fusion[k] <- "-1"
        }
      }
    }
  }
  all_net_result

}
