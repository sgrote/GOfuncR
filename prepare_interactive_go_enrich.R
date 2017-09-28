
# just some code to copy before going through go_enrich manually

setwd("/r1/people/steffi_grote/R_packages/GOfuncR_package")

library(GOfuncR)
library(data.table)
library(gtools)
library(mapplots)
load("GOfuncR/R/sysdata.rda")
source("GOfuncR/R/RcppExports.R")
source("GOfuncR/R/get_genes_from_regions.R")


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

# plot scores
fwer_threshold=0.26
go_ids=NULL
high_score_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_score_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
set.seed(123)
scores = c(sample(20:30, length(high_score_genes)), sample(5:15, length(low_score_genes)))
genes_wilcox = data.frame(genes_wilcox=c(high_score_genes, low_score_genes), scores)
res = go_enrich(genes_wilcox, test='wilcoxon')

# 2x2 contingency
test="contingency"
set.seed(123)
high_substi_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_substi_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2', 'C21orf59')
subs_syn = sample(45:55, length(c(high_substi_genes, low_substi_genes)), replace=T)
subs_non_syn = c(sample(15:25, length(high_substi_genes), replace=T), sample(0:10, length(low_substi_genes)))
vari_syn = sample(25:35, length(c(high_substi_genes, low_substi_genes)), replace=T)
vari_non_syn = c(sample(0:10, length(high_substi_genes), replace=T), sample(10:20, length(low_substi_genes)))
genes = data.frame(genes=c(high_substi_genes, low_substi_genes), vari_syn, vari_non_syn, subs_syn, subs_non_syn)
# plot odds ratio contingency
res = go_enrich(genes, test=test)
fwer_threshold = 0.7
go_ids = NULL


# binomial
set.seed(123)
high_human_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_human_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
human_counts = c(sample(20:30, length(high_human_genes)), sample(5:15, length(low_human_genes)))
chimp_counts = c(sample(5:15, length(high_human_genes)), sample(20:30, length(low_human_genes)))
genes = data.frame(gene=c(high_human_genes, low_human_genes), chimp_counts, human_counts)
# plot binom
res = go_enrich(genes, test='binomial', n_randsets=100)
fwer_threshold = 0.1
go_ids = NULL





