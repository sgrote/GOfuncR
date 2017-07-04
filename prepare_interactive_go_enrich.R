
# just some code to copy before going through go_enrich manually

setwd("/r1/people/steffi_grote/R_packages/FuncBlocks_package")

library(FuncBlocks)
library(data.table)
library(gtools)
load("FuncBlocks/R/sysdata.rda")
source("FuncBlocks/R/RcppExports.R")


# parameters

test="hyper"
n_randsets=100
gene_len=FALSE
circ_chrom=FALSE
ref_genome="grch37"
silent=FALSE
domains=NULL


gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
  'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
genes
