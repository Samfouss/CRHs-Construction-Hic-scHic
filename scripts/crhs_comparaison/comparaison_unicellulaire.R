load("rdata/all_rda_data/all_net_result_complex_3Mb_.rda")
load("rdata/all_rda_data/ideal_crhs_comparaison_res.rda")


lines = 0
for (bl in 2:16) {
  lines = length(all_net_result_complex_[[bl]]$more_info) + lines
}

# Function to process each line
process_line <- function(line) {
  # Step 1: Split line into parts separated by " - "
  parts <- unlist(strsplit(line, " - "))
  
  # Step 2: Extract numbers before '|' in each part
  cells <- sapply(parts, function(part) {
    return(as.numeric(gsub("[^0-9]", "", unlist(strsplit(part, "\\|"))[1])))
  })
  return(c(length(unique(cells)), length(cells)))
}

# Initialize matrix to store results
result_matrix <- matrix(NA, nrow = lines, ncol = 4)
colnames(result_matrix) <- c("Block", "Crh", "Numbre de cellules unique", "Number de cellule")

l = 1
for (bl in 2:16) {
  mat = all_net_result_complex_[[bl]]$more_info
  nb_crh = all_net_result_complex_[[bl]]$crhs
  for (crh in seq_along(nb_crh)) {
    result_matrix[l, 1] <- bl
    result_matrix[l, 2] <- crh
    result_matrix[l, 3] <- process_line(mat[crh])[1]
    result_matrix[l, 4] <- process_line(mat[crh])[2]
    l = l + 1
  }
}

# Display result
View(result_matrix)

table(result_matrix[, 1])
table(result_matrix[, 2])



# Analyse de la répetition des celluls et la sensibilité

sens_mat = crhs_comparaison_res$sensibility_mat
spec_mat = crhs_comparaison_res$specificity_mat

crhs_name = row.names(sens_mat)

rmin <- apply(sens_mat, 1, min, na.rm = TRUE)
rmean <- apply(sens_mat, 1, mean, na.rm = TRUE)
rmax <- apply(sens_mat, 1, max, na.rm = TRUE)

rep_cells_sens_distribution <- matrix(NA, nrow = length(crhs_name), ncol = 5)
colnames(rep_cells_sens_distribution) <- c("Numbre de cellules unique", "Number de cellule", "Min", "Moy", "Max")

for (crh in seq_along(crhs_name)) {
  text = crhs_name[crh]
  id = as.numeric(unlist(regmatches(text, gregexpr("[0-9]+", text))))
  crh_info = result_matrix[result_matrix[, 1]==id[1] & result_matrix[, 2]==id[2], ]
  # row_means <- rowMeans(df, na.rm = TRUE)
  rep_cells_sens_distribution[crh, 1] = crh_info[3]
  rep_cells_sens_distribution[crh, 2] = crh_info[4]
  rep_cells_sens_distribution[crh, 3] = rmin[crh]
  rep_cells_sens_distribution[crh, 4] = rmean[crh]
  rep_cells_sens_distribution[crh, 5] = rmax[crh]
  # apply(df, 1, max, na.rm = TRUE)
}

View(sens_mat)
View(rep_cells_sens_distribution)

summary(rep_cells_sens_distribution[rep_cells_sens_distribution[, 1]==1, 5])
summary(rep_cells_sens_distribution[rep_cells_sens_distribution[, 1]==2, 5])
summary(rep_cells_sens_distribution[rep_cells_sens_distribution[, 1]==3, 5])
summary(rep_cells_sens_distribution[rep_cells_sens_distribution[, 1]==4, 5])
summary(rep_cells_sens_distribution[rep_cells_sens_distribution[, 1]==5, 5])
summary(rep_cells_sens_distribution[rep_cells_sens_distribution[, 1]==6, 5])











