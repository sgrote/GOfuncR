
#test 'get_names', 'get_child_nodes', 'get_parent_nodes', 'get_anno_genes'

library(FuncBlocks)

##### get names
get_names(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
get_names('GO:0051082')
get_names(c('GO:123'))
get_names(c('bla','bla'))

##### get ids
get_ids('gaba')
get_ids('123')
get_ids(c('gaba','neuron'))
head(get_ids(''))
nrow(get_ids(''))

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
gos = c('GO:0072025','GO:0072221','GO:0072205','GO:0072235')
genes = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
anno_genes = get_anno_genes(go_ids=gos, genes=genes)
anno_genes
# get all genes annotated to these GOs
anno_all = get_anno_genes(go_ids=gos)
anno_all
# GO without annotation is not in output
get_anno_genes(go_ids=c('GO:007','GO:0072025'))
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
get_anno_genes(go_ids=c('GO:0072025','bla'))
get_anno_genes(go_ids=c('bla','GO:0072025'))
get_anno_genes(go_ids=c('GO:007'))


