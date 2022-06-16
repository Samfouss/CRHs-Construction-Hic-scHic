

find_min_by_rank <- function(element, min_rank){
  
  min_vec = c(min(element, na.rm = FALSE))
  
  if(min_rank > 1){
    for (i in seq_len(min_rank)) {
      print(element[min_vec != element])
      r_min = min(element[min_vec != element], na.rm = FALSE)
      min_vec = c(min_vec, r_min)
    } 
  }else{
    r_min = min_vec
  }
  
  r_min
}
