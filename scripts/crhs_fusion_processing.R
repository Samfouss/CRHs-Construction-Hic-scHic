library(ggplot2)
library(tidyverse)

load("rdata/all_net_result.rda")


result <- all_net_result$block1$resume_fusion%>%
  str_split("-", simplify = TRUE)

result <- gsub(",.*","",result)
result <- gsub(".* ","",result)

view(result)

result_df <- matrix(
  c(result[1, ][result[1, ]!= ""], str_c("crhs", rep(1, length(result[1, ][result[1, ]!= ""])))),
  nrow = length(result[1, ][result[1, ]!= ""]),
  ncol = 2,
  byrow = FALSE,
  dimnames = list(c(),c("CHR_NUM", "CRH_NAME"))
)

for (r in 2:nrow(result)) {
  
  result_df <- rbind(
    result_df,
    matrix(
      c(result[r, ][result[r, ]!= ""], str_c("crhs", rep(r, length(result[r, ][result[r, ]!= ""])))),
      nrow = length(result[r, ][result[r, ]!= ""]),
      ncol = 2,
      byrow = FALSE,
      dimnames = list(c(),c("CHR_NUM", "CRH_NAME"))
    )
  )
}

result_df <- data.frame(result_df)












