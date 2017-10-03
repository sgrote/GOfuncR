
# notepad on other GO-enrichment stuff on bioconductor



### topGO
library(topGO)
# example data
library(ALL)
data(ALL)
data(geneList) # named vector: p-values from differential expression named with genes

# apprently the user needs to provide the annotations, here included in the ALL-package:
affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)  # funktioniert nicht.. stop here


### CompGO
# can be used to download genome databases from UCSC entries
# interfaces DAVID which seems not to be maintained anymore

### EGSEA
# looks complicated at a first view

### geecc
# also looks complicated


### GOpro
# uses org.Hs.eg.db for annotations
# also to find groups of genes
# apparently only one domain can be run at a time

### GOSim
# uses topGO for the enrichment, also uses GO.db

### GOstats
# top 5%
# for Affymetrix arrays?
# one domain at a time
