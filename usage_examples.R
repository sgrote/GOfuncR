
set.seed(123)

###### standard go-enrichment with 14 candidate genes
library(FuncBlocks)

## input: a vector with '1' for candidate and optional '0' for background genes
## the names of the vector are the corresponding gene symbols	
gene_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
genes

# standard analysis with 1000 randomsets 
go_res = go_enrich(genes)
# the output is a list of 2 elements: 
# 1) the results from the anlysis (ordered by FWER for overrepresentation of candidate genes)
head(go_res[[1]])
by(go_res[[1]], go_res[[1]][,'ontology'], head)
# 2) the usable input genes 
go_res[[2]]



#### variations

# random gene selection dependent on gene-length (for FWER randomsets)
go_len = go_enrich(genes, gene_len=TRUE)

# different number of randomsets
go_less_ran = go_enrich(genes, n_randsets=100)

# define background
candi_gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
bg_gene_ids = c('FGR', 'NPHP1', 'DRD2', 'ABCC10', 'PTBP2', 'JPH4', 'SMARCC2', 'FN1', 'NODAL', 'CYP1A2', 'ACSS1', 'CDHR1', 'SLC25A36', 'LEPR', 'PRPS2', 'TNFAIP3', 'NKX3-1', 'LPAR2', 'PGAM2', 'GAPDHS')
genes = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
names(genes) = c(candi_gene_ids, bg_gene_ids)
genes
go_bg = go_enrich(genes, n_randsets=100)

# scores for genes and wilcoxon rank sum test
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = sample(1:30, length(gene_ids))
names(genes) = gene_ids
genes
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)


# regions as input
genes = c(1, rep(0,6))
names(genes) = c('8:82000000-83000000', '7:1300000-56800000', '7:74900000-148700000',
 '8:7400000-44300000', '8:47600000-146300000', '9:0-39200000', '9:69700000-140200000')
go_region = go_enrich(genes, n_randsets=100)
# genes located in candidate region
go_region[[2]][go_region[[2]]==1]

# background regions only on a circularized chromosome
go_region_circ = go_enrich(genes, n_randsets=100, circ_chrom=TRUE)





####### GO-graph functions (more to come)

#### (1) get gene annotations

## A) given a go_enrich result for GO-IDs with FWER<0.05 (FWER for overrepresentation / high_rank)
anno = get_anno_genes(go_res)
# using a different FWER-threshold
anno_willi = get_anno_genes(go_willi, fwer_threshold=0.2) 
# 'score' refers to candidate(1) or background(0) gene (or score in the wilcoxon test)
anno
anno_willi

# including background genes and using a different FWER-threshold
anno_bg = get_anno_genes(go_res, background=TRUE)
head(anno_bg)

# also works with genomic regions as input
anno_region = get_anno_genes(go_region, fwer_threshold=0.1)


## B) given GOs and 'optionally' genes directly
gos = c('GO:0072025','GO:0072221','GO:0072205','GO:0072235')
genes = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
anno_genes = get_anno_genes(go_ids=gos, genes=genes)
anno_genes

# get all genes annotated to these GOs
anno_all = get_anno_genes(go_ids=gos)
anno_all
# does not work with regions 
# (maybe 'yet', but its easy to get genes from regions and use those)
# (go_region[[2]] e.g. has all genes contained in the input regions for that go_enrich analysis)



#### (2) get name 
get_names(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))

### (3) GO->children
children = get_child_nodes(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
head(children)

### (4) GO->parents
parents = get_parent_nodes(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
parents

### (5) go_enrich-output and (fwer_threshold or go_ids)-> plot odds-ratios (hyper)
plot_odds_ratio(go_res, fwer_threshold=0.02)
plot_odds_ratio(go_bg,fwer_threshold=0.8)
plot_odds_ratio(go_res, go_ids=c('GO:0072025','GO:0072221','GO:0072235', 'GO:0044765'))
plot_odds_ratio(go_bg, go_ids=c('GO:0005634','GO:0004945','0.05309471','GO:0008289','GO:0005737','GO:0071495'))




#### future functions:
# (6) go_enrich-output and (fwer_threshold or go_ids)-> plot score-distribution (wilcoxon)








