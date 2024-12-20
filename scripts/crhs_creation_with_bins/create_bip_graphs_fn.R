

################# function 2 ##################
create_bip_graphs <- function(dt_structures, promoters_vec, rep_num, block_vec, min_dist = 3, cell_number, resolution=""){
  
  # Initilisation des des régions
  if(resolution == "1Mb"){
    resolution_part = c(1, 375, 750, 1125, 1500, 1875, 2248)
  }else if(resolution == "2Mb"){
    resolution_part = c(1, 749, 1499, 2248)
  }else if(resolution == "3Mb"){
    resolution_part = c(1, 1124, 2248)
  }else{
    resolution_part = c(1, 2248)
  }
  
  dt_structures <- dt_structures%>%
    mutate(
      is_promoter = if_else(row_number() %in% promoters_vec , 1, 0)
    )
  
  par_default <- par(bty = 'n')
  
  res = list(
    "block_1" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_2" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_3" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_4" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_5" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_6" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_7" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_8" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_9" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_10" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_11" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_12" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_13" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_14" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_15" = list(
      "more_info" = list(),
      "crhs" = list()
    ),
    "block_16" = list(
      "more_info" = list(),
      "crhs" = list()
    )
  )
  
  resolution_ = 0

  for (part in seq_len(length(resolution_part)-1)) {
    
    dt_structures_ = dt_structures[resolution_part[part]:resolution_part[part+1], ]
    block <- unique(pull(dt_structures_[, 4]))
    # print(block)
    # structure = data.frame()
    resolution_ = resolution_ + 1
    
    for(i in 1:length(block)) {
      
      b = block[i]
      
      three_dim_data <- dt_structures_[dt_structures_$X4==b, 1:3]
      block_promoters <- which(
        dt_structures_[dt_structures_$X4==b, "is_promoter"]$is_promoter==1
      )
      
      # S'il n y a pas de promoteurs dans le block
      if(length(block_promoters)!=0){
        
        mat_row_name = pull(dt_structures_[dt_structures_$X4==b, 5])
        #three_dim_data = column_to_rownames(three_dim_data, var = "ID")
        
        compute_dist <- as.matrix(
          dist(
            x = three_dim_data, 
            method = "euclidean", 
            diag = TRUE
          )
        )
        
        compute_dist_bin <- sapply(
          compute_dist, 
          # Si on est pas au block 1, on peut passer avec les contacts
          FUN = function(item) ifelse(item<=min_dist && b != 1, 1, 0)
        )
        
        
        compute_dist_bin <- matrix(
          compute_dist_bin, 
          nrow = nrow(compute_dist), 
          ncol = ncol(compute_dist), 
          byrow = TRUE,
          dimnames = list(mat_row_name, mat_row_name)
        )
        
        diag(compute_dist_bin) <- 0
        
        compute_dist_bin <- compute_dist_bin[block_promoters, -block_promoters, drop=FALSE]
        
        net_bip <- graph_from_incidence_matrix(
          compute_dist_bin
        )
        
        if(!is_simple(net_bip)){
          net_bip <- simplify(net_bip, remove.multiple = FALSE, remove.loops = TRUE)
        }
        net_components_bip <- components(net_bip, mode = c("weak", "strong"))
        # net_components <- components(net, mode = c("weak", "strong"))
        #biggest_cluster <- which.max(net_components$csize)
        #net_comp2 <- clusters(net, mode="weak")
        plot_main_cluster <- which(net_components_bip$csize>1)
        vert_ids <- V(net_bip)[net_components_bip$membership %in% plot_main_cluster]
        net_to_plot <- induced_subgraph(net_bip, vert_ids)
        V(net_to_plot)$color <- V(net_to_plot)$type
        V(net_to_plot)$color=gsub("FALSE","red",V(net_to_plot)$color)
        V(net_to_plot)$color=gsub("TRUE","lightblue",V(net_to_plot)$color)
        
        # A ce niveau, on prend chaque matrice d'incidence des différents CRHs qu'on ajoute à une liste 
        crhs = list()
        l1 = length(res[[b]]$crhs)
        l2 = length(res[[b]]$more_info)
        l = 1
        for (i in (1:net_components_bip$no)[net_components_bip$csize>1]){
          
          members <- net_components_bip$membership
          mat_bin = compute_dist_bin[
            rownames(compute_dist_bin) %in% names(members[members==i]),
            colnames(compute_dist_bin) %in% names(members[members==i]),
            drop=FALSE
          ]
          
          res[[b]]$crhs[[l1+l]] = list(
            "name" = str_c("( ", cell_number, "|", rep_num, ", ", i, ", ", resolution_, " )"),
            "mat_incidence" = mat_bin
          )
          l = l + 1
        }
        
        # Il y a des blocks qui ont des promoteurs par contre aucune bille n'est connetées à ce promoteur
        
        if(length(res[[b]]$crhs)>l1){
          names(res[[b]]$crhs) = str_c("crh", 1:length(res[[b]]$crhs))
          
          res[[b]]$more_info[[l2+1]] = list(
            "dist_bin" = compute_dist_bin,
            "net_to_plot" = net_to_plot,
            "membership_bip" = net_components_bip$membership,
            "csize_bip" = net_components_bip$csize,
            "no_bip" = net_components_bip$no
          )
          
        }
      }
    }
  }
  # print(length(result_name))
  names(res) <- str_c("block", 1:16)
  res
}




# cell = 10
# chr = 2
# data = as_tibble(
#   all_paired_structure%>%
#     filter(
#       paire == str_c(chr, sprintf("%03d", cell))
#     )%>%
#     select(-c(ends_with("_c"), "paire"))%>%
#     mutate(
#       ID = paste0("B", sprintf("%02d", X4), sprintf("%04d", 1:n()))
#     )
# )
# 
# 
# 
# all_net_result <- create_bip_graphs(
#   data,
#   promoters_ids,
#   chr,
#   1:16,
#   3,
#   sprintf("%03d", cell)
# )
# 
# ncrhs = 0
# for (bl in 1:16) {
#   print(length(all_net_result[[bl]]$crhs))
#   ncrhs = ncrhs + length(all_net_result[[bl]]$crhs)
# }
# ncrhs


