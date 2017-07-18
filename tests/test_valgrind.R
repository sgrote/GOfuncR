
## normal examples to run with valgrind

require(FuncBlocks)

##### hypergeometric standard parameters
set.seed(123)
gene_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
res = go_enrich(genes, n_randsets=100)
## gene len
res_len = go_enrich(genes, n_randsets=100, gene_len=T)

##### genomic regions
set.seed(123)
genes = c(1,1, rep(0,6))
names(genes) = c('8:82000000-83000000', '3:76500000-90500000', '7:1300000-56800000', '7:74900000-148700000',
 '8:7400000-44300000', '8:47600000-146300000', '9:0-39200000', '9:69700000-140200000')
go_region = go_enrich(genes, n_randsets=100)

##### wilcoxon
set.seed(123)
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = sample(1:30, length(gene_ids))
names(genes) = gene_ids
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)

##### binomial
set.seed(123)
high_human_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_human_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
human_counts = c(sample(20:30, length(high_human_genes)), sample(5:15, length(low_human_genes)))
chimp_counts = c(sample(5:15, length(high_human_genes)), sample(20:30, length(low_human_genes)))
genes = data.frame(gene=c(high_human_genes, low_human_genes), chimp_counts, human_counts)
go_binom = go_enrich(genes, test='binomial', n_randsets=100)

