
select_crh_candidats <- function(overlap_crh_object, overlap_edge_object) {
  
  results = list()
  
  for (b in 1:16){
    
    crh_match = overlap_crh_object[[b]]$chev_crh_comp
    edge_match = overlap_edge_object[[b]]$chev_edge_comp
    
    crh_edge_match = matrix(NA, 2, ncol = min(ncol(crh_match), ncol(edge_match)))
    
    for (m1 in seq_len(ncol(crh_match))) {
      for (m2 in seq_len(ncol(edge_match))) {
        if(all(crh_match[, m1]==edge_match[, m2])){
          crh_edge_match[, sum(!is.na(crh_edge_match[1, ])) + 1] <- crh_match[, m1]
        }
      }
    }
    
    crh_edge_match <- matrix(crh_edge_match[!is.na(crh_edge_match)], 2, ncol = length(crh_edge_match[!is.na(crh_edge_match)])/2)
    
    results[[length(results)+1]] <- list(
      "crh_edge_match" = crh_edge_match
    )
    
  }
  
  names(results) <- str_c("edge_overlap", 1:16)
  results
  
}


