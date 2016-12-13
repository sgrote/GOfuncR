
set.seed(123)

###### standard go-enrichment with 13 candidate genes
library(FuncBlocks)

## input: a vector with '1' for candidate and optional '0' for background genes
## the names of the vector are the corresponding gene symbols	
genes = rep(1, 13)
names(genes) = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')

# standard analysis with 1000 randomsets 
go_res = go_enrich(genes)
# the output is a list of 2 elements: the results from the anlysis and the input genes 
head(go_res[[1]]) # FWER_overrep is the most interesting
go_res[[2]] 



#### variations

# random gene selection dependent on gene-length (for FWER randomsets)
go_len = go_enrich(genes, gene_len=TRUE)

# different number of randomsets
go_less_ran = go_enrich(genes, n_randsets=100)

# scores for genes and wilcoxon rank sum test
genes = sample(1:30, 13)
names(genes) = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)


# regions as input
genes = c(1, rep(0,6))
names(genes) = c('8:82000000-83000000', '7:1300000-56800000', '7:74900000-148700000',
 '8:7400000-44300000', '8:47600000-146300000', '9:0-39200000', '9:69700000-140200000')
go_region = go_enrich(genes, n_randsets=100)

# background regions only on a circularized chromosome
go_region_circ = go_enrich(genes, n_randsets=100, circ_chrom=TRUE)



####### GO-graph functions (more to come)

#### (1) get gene annotations

## A) given a go_enrich result for GO-IDs with FWER<0.05 (FWER for overrepresentation / high_rank)
anno = get_anno_genes(go_res)
# using a different FWER-threshold
anno_willi = get_anno_genes(go_willi, fwer_threshold=0.2) 
# 'score' refers to candidate(1) or background(0) gene (or score in the wilcoxon test)
head(anno)
head(anno_willi)

# including background genes and using a different FWER-threshold
anno_bg = get_anno_genes(go_res, background=TRUE)
head(anno_bg)

# also works with genomic regions as input
anno_region = get_anno_genes(go_region, fwer_threshold=0.1)


## B) given GOs and "optionally" genes directly
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
get_names(c("GO:0051082", "GO:123", "GO:0042254", "GO:0000109"))

### (3) GO->children
children = get_child_nodes(c("GO:0051082", "GO:123", "GO:0042254", "GO:0000109"))
head(children)

### (4) GO->parents
parents = get_parent_nodes(c("GO:0051082", "GO:123", "GO:0042254", "GO:0000109"))
parents


#### future functions:
# (5) go_enrich-output,GO (or GO, genes+scores)-> plot odds-ratios (hyper) or score-distribution (wilcoxon)







