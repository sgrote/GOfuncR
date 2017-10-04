# install
# system('R CMD INSTALL /r1/people/steffi_grote/R_packages/GOfuncR_[VERSION].tar.gz')

# load
library(GOfuncR)

# random number seed
set.seed(123)


###### standard go-enrichment with 14 candidate genes

## input: a data-frame with two columns: (1) gene-identifiers and (2) scores: '1' for candidate and optional '0' for background genes, for the default hypergeometric test
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
gene_scores = rep(1, length(gene_ids))
genes = data.frame(gene_ids, gene_scores)
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
go_len = go_enrich(genes, n_randsets=100, gene_len=TRUE)

# only look at GO-domains "biological process" and "cellular component"
# (to save computation time)
go_2dom = go_enrich(genes, n_randsets=100, domains=c('biological_process','cellular_component'))

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
set.seed(123)
high_score_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_score_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = c(sample(20:30, length(high_score_genes)), sample(5:15, length(low_score_genes)))
names(genes) = c(high_score_genes, low_score_genes)
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)

# two counts for genes and binomial test
high_human_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_human_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
human_counts = c(sample(20:30, length(high_human_genes)), sample(5:15, length(low_human_genes)))
chimp_counts = c(sample(5:15, length(high_human_genes)), sample(20:30, length(low_human_genes)))
genes = data.frame(gene=c(high_human_genes, low_human_genes), human_counts, chimp_counts)
go_binom = go_enrich(genes, test='binomial', n_randsets=100)
head(go_binom[[1]])
go_binom[[2]]

# four counts per gene and 2x2-contingency-table test
# i.e. McDonald-Kreitman test where high rate of fixed non-synomymous changes / fixed synonymous changes between species compared to the rate of non-synonymous / synonymous variation inside one species could indicate positive selection
high_substi_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_substi_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2', 'C21orf59')
subs_non_syn = c(sample(15:25, length(high_substi_genes), replace=T), sample(0:10, length(low_substi_genes)))
subs_syn = sample(45:55, length(c(high_substi_genes, low_substi_genes)), replace=T)
vari_non_syn = c(sample(0:10, length(high_substi_genes), replace=T), sample(10:20, length(low_substi_genes)))
vari_syn = sample(25:35, length(c(high_substi_genes, low_substi_genes)), replace=T)
genes = data.frame(genes=c(high_substi_genes, low_substi_genes), subs_non_syn, subs_syn, vari_non_syn, vari_syn)
conti_res = go_enrich(genes, test='contingency', n_randset=100)
head(conti_res[[1]])

# regions as input
genes = c(1, rep(0,6))
names(genes) = c('8:81000000-83000000', '7:1300000-56800000', '7:74900000-148700000',
 '8:7400000-44300000', '8:47600000-146300000', '9:0-39200000', '9:69700000-140200000')
genes
go_region = go_enrich(genes, n_randsets=100)
# genes located in candidate region
go_region[[2]][go_region[[2]][,2]==1,]

# use hg20 or mouse gene coordinates to find genes located in the regions
go_region_hg20 = go_enrich(genes, n_randsets=100, ref_genome='grch38')
go_region_mus = go_enrich(genes, n_randsets=100, ref_genome='grcm38')

# background regions only on a circularized chromosome
go_region_circ = go_enrich(genes, n_randsets=100, circ_chrom=TRUE)


## since version 1.2.4 parallel processing is possible
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




####### GO-graph functions


#### (1) GO-ID -> genes  
# (NEW: this function does not support go_enrich()-result as input anymore (simpler))

# find all genes that are annotated to GO:0000109
# ("nucleotide-excision repair complex")
get_anno_genes(go_ids='GO:0000109', ref_genome='grch38')

# find out wich genes from a set of genes
# are annotated to some GO-categories
genes = c('AGTR1', 'ANO1', 'CALB1', 'GYG1', 'PAX2')
gos = c('GO:0001558', 'GO:0005536', 'GO:0072205', 'GO:0006821')
anno_genes = get_anno_genes(go_ids=gos, genes=genes)
# and add the names and domains of the GO-categories
cbind(anno_genes, get_names(anno_genes$go_id)[,2:3])

# find all mouse-gene annotations to two GO-categories 
gos = c('GO:0072205', 'GO:0000109')
get_anno_genes(go_ids=gos, ref_genome='grcm38')

# extract categories with a FWER<0.05 from an enrichment analysis
# and find the candidate genes annotated to those top-hits
# (this works the same way for an analysis with genomic regions)
stats = go_res[[1]]
input_genes = go_res[[2]]
candidate_genes = input_genes[input_genes[,2]==1, 1]
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

#### (3) GO-name -> GO-ID  (all partial perfect matches of the input string, not case sensitive) 
get_ids(c('ribosome'))
head(get_ids(c('gaba')))

#### (4) GO-ID -> children
children = get_child_nodes(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
head(children)

#### (5) GO-ID -> parents
parents = get_parent_nodes(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
parents

#### (6) plot_anno_score






