
############ Programme pour retirer les bins aléatoirement du bloc 1 #######

rep_number = 50

fs <- list.files(path = 'rdata/', pattern = '.txt')
for (f in fs) {
  
  path<-paste("rdata/",f,sep="")
  structure <- read.table(path)
  num_to_remove = sample(seq_len(nrow(structure)), 2, replace = FALSE)
  
  structure <- structure[-num_to_remove,]
  # do all your processing stuff to df
  write.table(
    x = structure, 
    file = f, 
    fileEncoding = "UTF-8",
    col.names=FALSE,
    row.names=FALSE
  )

}

