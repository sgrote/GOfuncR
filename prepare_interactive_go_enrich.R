
# just some code to copy before going through go_enrich manually

setwd("/r1/people/steffi_grote/R_packages/FuncBlocks_package")

library(FuncBlocks)
library(data.table)
library(gtools)
load("FuncBlocks/R/sysdata.rda")
source("FuncBlocks/R/RcppExports.R")


# parameters
n_randsets=100
gene_len=FALSE
circ_chrom=FALSE
ref_genome="grch37"
silent=FALSE
domains=NULL

# hyper
test="hyper"
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
  'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
genes

# 2x2 contingency
test="contingency"
gene_ids = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE', 'CACNG2')
subs_syn = sample(45:55, length(gene_ids), replace=T)
subs_non_syn = sample(15:25, length(gene_ids), replace=T)
vari_syn = sample(25:35, length(gene_ids), replace=T)
vari_non_syn = sample(0:10, length(gene_ids), replace=T)
genes = data.frame(gene_ids, subs_syn, subs_non_syn, vari_syn, vari_non_syn)
