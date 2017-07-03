# install
system('R CMD INSTALL /r1/people/steffi_grote/R_packages/FuncBlocks_1.2.6.tar.gz')

# load
library(FuncBlocks)

# random number seed
set.seed(123)


###### standard go-enrichment with 14 candidate genes

## input: a vector with '1' for candidate and optional '0' for background genes
## the names of the vector are the corresponding gene symbols	
gene_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
genes

# standard analysis with 1000 randomsets 
go_res = go_enrich(genes)
# the output is a list of 3 elements: 
# 1) the results from the anlysis (ordered by FWER for overrepresentation of candidate genes)
head(go_res[[1]])
by(go_res[[1]], go_res[[1]][,'ontology'], head)
# 2) all valid input genes
go_res[[2]]
# 3) reference genome used (default=grch37)
go_res[[3]]


#### variations

# different number of randomsets (NOTE: also set to 100 in the following examples to lower computation time)
go_less_ran = go_enrich(genes, n_randsets=100)

# random gene selection dependent on gene-length (for FWER randomsets)
go_len = go_enrich(genes, gene_len=TRUE)

# grch38 reference genome GO-annotations (default=grch37)
go_hg20 = go_enrich(genes, ref_genome='grch38', n_randsets=100)

# mouse reference genome
gene_ids = c('Arsi', 'Mapk4', 'Papola', 'Tfrc', 'Bak1', 'Fopnl', 'Mus81', 'Opa3', 'Npcd')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
go_mouse = go_enrich(genes, ref_genome='grcm38', n_randsets=100)

# define background
set.seed(123)
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
names(genes) = c('8:81000000-83000000', '7:1300000-56800000', '7:74900000-148700000',
 '8:7400000-44300000', '8:47600000-146300000', '9:0-39200000', '9:69700000-140200000')
genes
go_region = go_enrich(genes, n_randsets=100)
# genes located in candidate region
go_region[[2]][go_region[[2]]==1]

# use hg20 or mouse gene coordinates to find genes located in the regions
go_region_hg20 = go_enrich(genes, n_randsets=100, ref_genome='grch38')
go_region_mus = go_enrich(genes, n_randsets=100, ref_genome='grcm38')

# background regions only on a circularized chromosome
go_region_circ = go_enrich(genes, n_randsets=100, circ_chrom=TRUE)


## NEW: since version 1.2.4 parallel processing is possible
library("parallel")
# create a list with 3 different input 'genes' vectors
candi1_gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'NPHP1', 'DRD2', 'FN1', 'NODAL')
candi2_gene_ids = c('C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'ABCC10', 'PTBP2', 'CYP1A2', 'ACSS1')
candi3_gene_ids = c('BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2', 'JPH4', 'SMARCC2', 'CDHR1', 'SLC25A36')
bg_gene_ids = c('FGR', 'LEPR', 'PRPS2', 'TNFAIP3', 'NKX3-1', 'LPAR2', 'PGAM2', 'GAPDHS')
input=list()
input[["genes1"]] = c(rep(1,length(candi1_gene_ids)), rep(0,length(bg_gene_ids)))
names(input[["genes1"]]) = c(candi1_gene_ids, bg_gene_ids)
input[["genes2"]] = c(rep(1,length(candi2_gene_ids)), rep(0,length(bg_gene_ids)))
names(input[["genes2"]]) = c(candi2_gene_ids, bg_gene_ids)
input[["genes3"]] = c(rep(1,length(candi3_gene_ids)), rep(0,length(bg_gene_ids)))
names(input[["genes3"]]) = c(candi3_gene_ids, bg_gene_ids)
# run go_enrich for all 3 input-vectors in parallel
parares = mclapply(1:3, function(x){
	set.seed(123)
	go_enrich(input[[x]], n_randset=50)
})




####### GO-graph functions (more to come)


#### (1) GO-ID -> genes  
# (NEW: since version 1.2.6 this function does not support go_enrich()-result as input anymore (simpler))

# find all genes that are annotated to GO:0000109
# ("nucleotide-excision repair complex")
get_anno_genes(go_ids='GO:0000109', ref_genome='grch38')

# find out wich genes from a set of genes
# are annotated to some GO-categories
genes = c('AGTR1', 'ANO1', 'CALB1', 'GYG1', 'PAX2')
gos = c('GO:0001558', 'GO:0005536', 'GO:0072205', 'GO:0006821')
anno_genes = get_anno_genes(go_ids=gos, genes=genes)
# and add the names and domains of the GO-categories
cbind(anno_genes ,get_names(anno_genes$go_id)[,2:3])

# find all mouse-gene annotations to two GO-categories 
gos = c('GO:0072205', 'GO:0000109')
get_anno_genes(go_ids=gos, ref_genome='grcm38')

# extract categories with a FWER<0.05 from an enrichment analysis
# and find the candidate genes annotated to those top-hits
stats = go_res[[1]]
candidate_genes = names(go_res[[2]][go_res[[2]]==1])
ref_genome = go_res[[3]]
top_hits = stats[stats$FWER_overrep < 0.05, 'node_id']
anno_top = get_anno_genes(go_ids=top_hits, ref_genome=ref_genome, genes=candidate_genes)
# and add the names and domains of the GO-categories
cbind(anno_top ,get_names(anno_top$go_id)[,2:3])

# extract categories with a FWER<0.05 from an enrichment analysis
# and find the candidate genes annotated to those top-hits
# (this works the same way for an analysis with genomic regions)
stats = go_res[[1]]
candidate_genes = names(go_res[[2]][go_res[[2]]==1])
ref_genome = go_res[[3]]
top_hits = stats[stats$FWER_overrep < 0.05, 'node_id']
anno_top = get_anno_genes(go_ids=top_hits, ref_genome=ref_genome, genes=candidate_genes)
# and add the names and domains of the GO-categories
cbind(anno_top ,get_names(anno_top$go_id)[,2:3])


#### (2) gene -> GO-IDs
get_anno_categories(c('NCAPG', 'APOL4'))
get_anno_categories(c('Mus81', 'Papola'), ref_genome='grcm38')

#### (3) GO-ID -> GO-name 
get_names(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))

#### (4) GO-ID -> children
children = get_child_nodes(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
head(children)

#### (5) GO-ID -> parents
parents = get_parent_nodes(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
parents

#### (6) go_enrich-output and (fwer_threshold or go_ids) -> plot odds-ratios (hyper)
plot_odds_ratio(go_res, fwer_threshold=0.02)
plot_odds_ratio(go_bg,fwer_threshold=0.8)
plot_odds_ratio(go_res, go_ids=c('GO:0072025','GO:0072221','GO:0072235', 'GO:0044765'))
plot_odds_ratio(go_bg, go_ids=c('GO:0005634','GO:0004945','0.05309471','GO:0008289','GO:0005737','GO:0071495'))




#### future functions:
# (7) go_enrich-output and (fwer_threshold or go_ids)-> plot score-distribution (wilcoxon)








