## Installation de bioconductor si c'est pas le cas sur la machine
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

## Installation des packages dont on aura probablement besoin
BiocManager::install("IRanges")
BiocManager::install("rtracklayer") # Ce package permet d'importer les fichiers génomiques de différents formats
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("org.Mm.eg.db")
## Chargeent des librairy
library(rtracklayer)
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
symbol <- select(
  org.Mm.eg.db,
  keys = genes_muscu$gene_id, 
  keytype = "ENTREZID", 
  columns = "SYMBOL"
)
genes_muscu$geneSymbol <- symbol$SYMBOL
genesInRRCs = genes_muscu[genes_muscu$geneSymbol%in%musculusSox9Gene$Name]

promoters = promoters(genesInRRCs)

df.promoters = data.frame(promoters)
df.promoters = df.promoters[, c("seqnames", "start", "end", "geneSymbol")]
colnames(df.promoters) = c("chr", "startProm", "endProm", "Name")

enhancers.promoters = merge(musculusSox9Gene, df.promoters, by="Name")




