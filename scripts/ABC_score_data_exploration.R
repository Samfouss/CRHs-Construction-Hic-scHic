library(tidyverse)

endo_cell_data_02 <- read_tsv("rdata/SCORE_ABC_RESULTS/ENDO_CELL/EnhancerPredictions_02.txt", show_col_types = FALSE)
summary(as.numeric(table(endo_cell_data_02$TargetGene)))


astro_cell_data_05 <- read_tsv("rdata/SCORE_ABC_RESULTS/ASTRO_CELL/EnhancerPredictions_05.txt", show_col_types = FALSE)
summary(as.numeric(table(astro_cell_data_05$TargetGene)))

astro_cell_data_03 <- read_tsv("rdata/SCORE_ABC_RESULTS/ASTRO_CELL/EnhancerPredictions_03.txt", show_col_types = FALSE)
summary(as.numeric(table(astro_cell_data_03$TargetGene)))

astro_cell_data_02 <- read_tsv("rdata/SCORE_ABC_RESULTS/ASTRO_CELL/EnhancerPredictions_02.txt", show_col_types = FALSE)
summary(as.numeric(table(astro_cell_data_02$TargetGene)))

astro_cell_data_01 <- read_tsv("rdata/SCORE_ABC_RESULTS/ASTRO_CELL/EnhancerPredictions_01.txt", show_col_types = FALSE)
summary(as.numeric(table(astro_cell_data_01$TargetGene)))

astro_cell_data_005 <- read_tsv("rdata/SCORE_ABC_RESULTS/ASTRO_CELL/EnhancerPredictions_005.txt", show_col_types = FALSE)
summary(as.numeric(table(astro_cell_data_005$TargetGene)))



