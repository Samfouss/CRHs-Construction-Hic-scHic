
create_match_merge_crhs <- function(all_structures){
  
  # On crée ici le tout premier CRHs à partir de la première structure
  result <- create_bip_clust_graph(
    all_structures%>%
      filter(replica == 1)%>%
      select("X1", "X2", "X3", "X4"), 
    promoters_ids, 
    1, 
    1:16, 
    3
  )
  
  for (str in 2:500) {
    
    structure_selected <- all_structures%>%
      filter(replica == str)%>%
      select("X1", "X2", "X3", "X4")
    
    structure_selected$ID <- paste0("B", sprintf("%02d", structure_selected$X4), sprintf("%04d", 1:nrow(structure_selected)))
    
    # Take 1 structure create a CRHs with
    # On peut poser une condition à la première itération pour ne pas faire la fusion, on garde les résultats à une liste "final_result"
    # Lors de la deuxième itération, on prend la structure suivante, on creé des CRHs avec puis on fusionne avec les CRHs de "final_result", le résultat sera restocké dans la liste "final_result" et progressivement
    
    net_bip_for_sected_structure <- create_bip_clust_graph(
      structure_selected, 
      promoters_ids, 
      str, 
      1:16, 
      3
    )
    
    overlap_edge <- edge_identity_overlap(
      result, 
      net_bip_for_sected_structure
    )
    
    result <- fusion_comp(result, net_bip_for_sected_structure, overlap_edge)
    
  }
  
  
}


