
# just some code to copy before going through go_enrich manually

setwd("/r1/people/steffi_grote/R_packages/FuncBlocks_package")

library(FuncBlocks)
library(data.table)
library(gtools)
load("FuncBlocks/R/sysdata.rda")
source("FuncBlocks/R/RcppExports.R")
source("FuncBlocks/R/get_genes_from_regions.R")


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

# genomic regions
test="hyper"
genes = c(1,1, rep(0,6))
names(genes) = c('8:81000000-83000000', '3:76500000-90500000', '7:1300000-56800000', '7:74900000-148700000', '8:7400000-44300000', '8:47600000-146300000', '9:0-39200000', '9:69700000-140200000')


# wilcoxon
test="wilcoxon"
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = sample(1:30, length(gene_ids))
names(genes) = gene_ids

# 2x2 contingency
test="contingency"
set.seed(123)
gene_ids = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE', 'CACNG2')
subs_syn = sample(45:55, length(gene_ids), replace=T)
subs_non_syn = sample(15:25, length(gene_ids), replace=T)
vari_syn = sample(25:35, length(gene_ids), replace=T)
vari_non_syn = sample(0:10, length(gene_ids), replace=T)
genes = data.frame(gene_ids, subs_syn, subs_non_syn, vari_syn, vari_non_syn)
# plot odds ratio contingency
res = go_enrich(genes, test=test)
fwer_threshold=0.26
go_ids=NULL


# plot scores
fwer_threshold=0.26
go_ids=NULL
high_score_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_score_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
set.seed(123)
scores = c(sample(20:30, length(high_score_genes)), sample(5:15, length(low_score_genes)))
genes_wilcox = data.frame(genes_wilcox=c(high_score_genes, low_score_genes), scores)
res = go_enrich(genes_wilcox, test='wilcoxon')










