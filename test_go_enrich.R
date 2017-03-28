# try out functionalities, error messages and warnings

# remove package from R session
detach('package:FuncBlocks', unload = TRUE)
#library.dynam.unload("FuncBlocks", system.file(package = "FuncBlocks"))

##############################

set.seed(123)
library(FuncBlocks)
setwd('/r1/people/steffi_grote/R_packages/FuncBlocks_package')


##### standard parameters
gene_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
res = go_enrich(genes)
head(res[[1]])
res[[2]] # should not contain QUATSCH1
### corner cases
# one gene
res = go_enrich(genes[1], n_randsets=100)
head(res[[1]])
# one gene without annotation
res = go_enrich(genes[2])
### erroneous input
# not 1/0-input
genes2 = genes
genes2[3] = 2
res = go_enrich(genes2)
# no genes names
res = go_enrich(c(1,0,0,1,0))
# no gene scores
res = go_enrich(gene_ids)


##### standard parameters - background defined
candi_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4')
bg_ids = c('C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = c(rep(1,length(candi_ids)), rep(0,length(bg_ids)))
names(genes) = c(candi_ids, bg_ids)
res = go_enrich(genes, n_randsets=100)
head(res[[1]])
res[[2]]
### corner cases
# one background gene
res = go_enrich(genes[1:(length(candi_ids)+1)], n_randsets=100)
head(res[[1]])
res[[2]]
# one candidate gene
res = go_enrich(genes[c(1,(length(candi_ids)+1):length(genes))], n_randsets=100)
head(res[[1]])
res[[2]]
### erroneous input
# same gene as candidate and background
names(genes)[1:3] = bg_ids[1:3]
res = go_enrich(genes)
# only background defined
res = go_enrich(genes[(length(candi_ids)+1):length(genes)])


##### wilcoxon
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = sample(1:30, length(gene_ids))
names(genes) = gene_ids
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)
head(go_willi[[1]])
go_willi[[2]]
### corner cases
# negative values
genes_rev=-genes
go_willi_rev = go_enrich(genes_rev, test='wilcoxon', n_randsets=100)
head(go_willi_rev[[1]])
forward_p_low = go_willi[[1]][match(go_willi_rev[[1]][,2], go_willi[[1]][,2]),'raw_p_low_rank']
all.equal(go_willi_rev[[1]][,'raw_p_high_rank'], forward_p_low)
# only two genes
go_willi = go_enrich(genes[1:2], test='wilcoxon', n_randsets=100)
head(go_willi[[1]])
go_willi[[2]]
# only one gene
go_willi = go_enrich(genes[3], test='wilcoxon', n_randsets=100)
# only two scores
genes = sample(1:2, length(gene_ids), replace=TRUE)
names(genes) = gene_ids
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)
head(go_willi[[1]])
# only one score - works, all p and FWER are 1
genes[genes==1] = 2
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)
head(go_willi[[1]])
# floating point input
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1')
genes = seq(1.1, 1.7, by=0.1)
names(genes) = gene_ids
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)
head(go_willi[[1]])


##### n_randsets
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
### corner cases
go_ran = go_enrich(genes, n_randsets=0)
go_ran = go_enrich(genes, n_randsets=1)
# float
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=10.5)
### erroneous input
go_enrich(genes, n_randsets='ahh')
go_enrich(genes, n_randsets=-3)


##### gene_len
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
res_len = go_enrich(genes, gene_len=TRUE, n_randset=100)
head(res_len[[1]])
### corner cases
# no input gene has coordinates
no_coord_ids = c('HMGA1P6', 'RNY3P4', 'LINC00362', 'RNU6-58P', 'TATDN2P3', 'LINC00363')
names(genes) = no_coord_ids
res_len = go_enrich(genes, gene_len=TRUE, n_randset=100)
# no background has coordinates
genes_bg = rep(c(1,0),each=6)
names(genes_bg) = c(gene_ids, no_coord_ids)
res_len = go_enrich(genes, gene_len=TRUE, n_randset=100)
### errorneous input
# other than T/F
res_len = go_enrich(genes, gene_len='bla', n_randset=100)
res_len = go_enrich(genes, gene_len=1:3, n_randset=100)
# wilcoxon
res_len = go_enrich(genes, test='wilcoxon', gene_len=TRUE, n_randset=100)


##### genomic regions
genes = c(1,1, rep(0,6))
names(genes) = c('8:82000000-83000000', '3:76500000-90500000', '7:1300000-56800000', '7:74900000-148700000',
 '8:7400000-44300000', '8:47600000-146300000', '9:0-39200000', '9:69700000-140200000')
go_region = go_enrich(genes, n_randsets=100)
head(go_region[[1]])
head(go_region[[2]])
### corner cases
# background too tight for random placement
tight = c(rep(1,4),0)
names(tight) = c('1:104000000-114900000', '3:76500000-90500000', '7:113600000-124700000', '8:54500000-65400000', '5:0-4700000')
go_enrich(tight, n_randsets=100)
### erroneous input
# start larger than end
reverse = tight
names(reverse)[3] = '2:15-10'
go_enrich(reverse)
# no background 
go_enrich(genes[1:2])
# no candi
go_enrich(genes[3:length(genes)])
# background too small
too_small = c(rep(1,4),rep(0,5)) 
names(too_small) = c('X:0-2','2:0-3','2:5-10','13:0-20',  '2:5-10','13:0-10','X:0-1','4:0-10','2:10-12') 
go_enrich(too_small)
# no genes in candidate
no_can = genes
names(no_can)[1:2] = c('1:10-20', '21:1000-3000')
go_enrich(no_can)
# no genes in background
no_bg = genes[c(1,3)]
names(no_bg)[2] = c('21:1-3000000')
go_no_bg = go_enrich(no_bg)


##### circ_chrom
go_circ = go_enrich(genes[c(1,5,6)], n_randsets=100, circ_chrom=TRUE)
head(go_circ[[1]])
# warning about unused chromosomes
go_circ_unused = go_enrich(genes[-2], n_randsets=100, circ_chrom=TRUE)
head(go_circ_unused[[1]])
# unused chromosomes not contained in returned genes
length(go_circ[[2]][go_circ[[2]]==0]) == length(go_circ_unused[[2]][go_circ_unused[[2]]==0])
length(go_circ_unused[[2]][go_circ_unused[[2]]==0]) < length(go_region[[2]][go_region[[2]]==0])

##### example with Ben's deserts
background = read.table('Ben_background_regions.bed')
candidate = read.table('Ben_neandertal_deserts.bed')
back_cand = rbind(candidate, background)
regions = c(rep(1,nrow(candidate)),rep(0,nrow(background)))
names(regions) = c(paste(back_cand[,1],':',back_cand[,2],'-',back_cand[,3],sep=''))
# normal blocks 
ben = go_enrich(regions)
# circ chrom
ben_circ = go_enrich(regions, circ_chrom=TRUE)
# single genes extracted from normal blocks
ben_genes = go_enrich(ben[[2]])
# single genes extracted from circ chrom
ben_genes_circ = go_enrich(ben_circ[[2]])
## p-values should be the same for single genes
all.equal(ben[[1]][,5],ben_genes[[1]][match(ben[[1]][,'node_id'],ben_genes[[1]][,'node_id']),5])
all.equal(ben_circ[[1]][,5],ben_genes_circ[[1]][match(ben_circ[[1]][,'node_id'],ben_genes_circ[[1]][,'node_id']),5])
## circ and non-circ should NOT match due to different background set
all.equal(ben_circ[[1]][,5],ben[[1]][match(ben_circ[[1]][,'node_id'],ben[[1]][,'node_id']),5])




