# Installation de bioconductor si c'est pas le cas sur la machine
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

## Installation des packages dont on aura probablement besoin
# BiocManager::install("IRanges")
# BiocManager::install("rtracklayer") # Ce package permet d'importer les fichiers génomiques de différents formats
# BiocManager::install("rtracklayer") # Ce package permet d'importer les fichiers génomiques de différents formats
# BiocManager::install("org.Mm.eg.db")
## Chargeent des librairy
library(rtracklayer)
library(IRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

## Importation des données
musculusSox9Gene <- import("rdata/musculusSox9GeneChr11-109–115.GFF3")

muscu_txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
keytypes(muscu_txdb)
keys(muscu_txdb)
keytypes(org.Mm.eg.db)
keys(org.Mm.eg.db)
genes_muscu <- genes(muscu_txdb)
symbol <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = genes_muscu$gene_id, 
  keytype = "ENTREZID", 
  columns = "SYMBOL"
)
genes_muscu$geneSymbol <- symbol$SYMBOL
genesInRRCs = genes_muscu[genes_muscu$geneSymbol%in%musculusSox9Gene$Name]

promoters = promoters(genesInRRCs)
# On filtre les promoters qui ne sont pas sur le chroosome 11
promotersChr11 = promoters[promoters@seqnames=="chr11"]


#################### Get promoters for scHic ######################
# Les matrices scHic contiennent par bac 4 billes.
# Pour une bille qui vaut 6000000/2248=2667, les bacs 4 billes vallent donc (6000000/2248)*4 = 10676
polymers_range <- seq(from = (109000000-1), to = 115000000, by = 10676)
polymers_query <- IRanges(
  polymers_range[1:(length(polymers_range))]+1,
  c(polymers_range[2:length(polymers_range)], 115000000),
  names=sprintf("%04d", 1:length(polymers_range))
)
# Récupération des IRanges sur le chromozome 11
polymers_subject <- promotersChr11@ranges
# Trouver les chevauchements
polymers_promoter <- findOverlaps(polymers_query, polymers_subject)
# Récupérer les IDs des billes sans les dupliquées
scHic_promoters_ids <- polymers_promoter@from[!duplicated(polymers_promoter@from)]

# Save promoters ID
save(scHic_promoters_ids, file = "rdata/scHic_promoters_ids.rda")




