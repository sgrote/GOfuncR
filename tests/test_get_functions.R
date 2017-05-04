
#test 'get_names', 'get_child_nodes', 'get_parent_nodes', 'get_anno_genes'

library(FuncBlocks)

##### get names
get_names(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
get_names('GO:0051082')
get_names(c('GO:123'))
get_names(c('bla','bla'))

##### get_child_nodes
get_child_nodes(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
get_child_nodes('GO:0051082')
get_child_nodes(c('GO:123'))
get_child_nodes(c('bla','bla'))

##### get_parent_nodes
get_parent_nodes(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
get_parent_nodes('GO:0003674') #root
get_parent_nodes(c('GO:123'))
get_parent_nodes(c('bla','bla'))


##### get_anno_genes

### A) given a go_enrich result for GO-IDs with FWER<0.05 (FWER for overrepresentation / high_rank)

set.seed(123)
# standard analysis with 1000 randomsets 
gene_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
go_res = go_enrich(genes)
# normal
anno = get_anno_genes(go_res)
anno_bg = get_anno_genes(go_res, background=TRUE)
anno_all = get_anno_genes(go_res, fwer_threshold=5) # its possible to get all annotations
# check that nr. of candidate genes matches output from go_enrich
anno_sum = tapply(anno$score, anno$go_id, sum)
all(anno_sum == go_res[[1]][match(names(anno_sum), go_res[[1]][,2]), 'n_candidate_real'])
# TODO: when output of total number of annotated genes is implemented in go_enrich -> test that as well
# no FWER below threshold
get_anno_genes(go_res, fwer_threshold=0.005) # 0.005=minimum fwer in res
# unused arguments: Warnings
get_anno_genes(go_res, go_ids='GO:0072025')
get_anno_genes(go_res, genes='CALB1')
get_anno_genes(go_res, go_ids='GO:0072025', genes='CALB1')

# background defined
candi_gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
bg_gene_ids = c('FGR', 'NPHP1', 'DRD2', 'ABCC10', 'PTBP2', 'JPH4', 'SMARCC2', 'FN1', 'NODAL', 'CYP1A2', 'ACSS1', 'CDHR1', 'SLC25A36', 'LEPR', 'PRPS2', 'TNFAIP3', 'NKX3-1', 'LPAR2', 'PGAM2', 'GAPDHS')
genes = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
names(genes) = c(candi_gene_ids, bg_gene_ids)
go_bg = go_enrich(genes, n_randsets=100)

# scores for genes and wilcoxon rank sum test
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = sample(1:30, length(gene_ids))
names(genes) = gene_ids
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)
anno_willi = get_anno_genes(go_willi, fwer_threshold=0.93)
# check that score matches input
all(anno_willi$score == genes[match(anno_willi[,2], names(genes))])

# regions as input
genes = c(1, rep(0,6))
names(genes) = c('8:82000000-83000000', '7:1300000-56800000', '7:74900000-148700000',
 '8:7400000-44300000', '8:47600000-146300000', '9:0-39200000', '9:69700000-140200000')
go_region = go_enrich(genes, n_randsets=100)
anno_region = get_anno_genes(go_region, fwer_threshold=0.1)
anno_region_bg = get_anno_genes(go_region, fwer_threshold=0.1, background=TRUE)


### B) given GOs and 'optionally' genes directly
gos = c('GO:0072025','GO:0072221','GO:0072205','GO:0072235')
genes = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
anno_genes = get_anno_genes(go_ids=gos, genes=genes)
anno_genes
# get all genes annotated to these GOs
anno_all = get_anno_genes(go_ids=gos)
anno_all
# no go_ids
get_anno_genes(genes=genes)
get_anno_genes(123)
get_anno_genes(c('bla','bla'))
# genes invalid -> no annotation
get_anno_genes(go_ids=gos, genes=123)
get_anno_genes(go_ids=gos, genes=c('bla','bla'))
# go_ids invalid
get_anno_genes(go_ids=123)
get_anno_genes(go_ids=c('bla','bla'))

