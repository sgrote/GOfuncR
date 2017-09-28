# try out functionalities, error messages and warnings

# remove package from R session
detach('package:GOfuncR', unload = TRUE)
#library.dynam.unload("GOfuncR", system.file(package = "GOfuncR"))

##############################

library(GOfuncR)
setwd('/r1/people/steffi_grote/R_packages/GOfuncR_package')


##### standard parameters
gene_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
scores = rep(1, length(gene_ids))
genes = data.frame(gene_ids, scores) 
res = go_enrich(genes)
head(res[[1]])
res[[2]] # should not contain QUATSCH1
### corner cases
# one gene
res = go_enrich(genes[1,], n_randsets=100)
head(res[[1]])
### erroneous input
# one gene without annotation
res = go_enrich(genes[2,])
# not 1/0-input
genes2 = genes
genes2[3,2] = 2
res = go_enrich(genes2)
# no genes names
res = go_enrich(genes[,2])
# no gene scores
res = go_enrich(genes[,1])



##### standard parameters - background defined
candi_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4')
bg_ids = c('C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
scores = c(rep(1,length(candi_ids)), rep(0,length(bg_ids)))
genes = data.frame(g=c(candi_ids, bg_ids), scores)
set.seed(123)
res = go_enrich(genes, n_randsets=100)
head(res[[1]])
res[[2]]
# vector_input
set.seed(123)
gene_vec = data.frame(a=names(genes), b=unname(genes))
res_vec = go_enrich(genes, n_randsets=100)
all.equal(res, res_vec)
# multiple assignment of same value - ok
res = go_enrich(genes[c(1,1:length(genes)),], n_randsets=10)
### corner cases
# one background gene (not annotated in all roots - only cellular component)
onebg = genes[c(1,(length(candi_ids)+1):length(genes)),]
res = go_enrich(onebg, n_randsets=100)
head(res[[1]])
by(res[[1]], res[[1]]$ontology, summary) # all p and FWER = 1 in mole and biol roots
res[[2]]
# one candidate gene (not annotated in all roots - only cellular component)
onecan = onebg
onecan[,2] = 1*(onecan[,2]==0)
res = go_enrich(onecan, n_randsets=100)
head(res[[1]])
by(res[[1]], res[[1]]$ontology, summary) # all p and FWER = 1 in mole and biol roots
res[[2]]
## candidate AND background lack annotations in a specific root?
# no candidate in mol/biol; no background in mol
ungenes = data.frame(a=c("C21orf59","MTUS1"), b=c(1,0))  ## TODO: this should be skipped in next version
he = go_enrich(ungenes, n_randset=50) ## vorher: Error in evalq(sys.calls(), <environment>
by(he[[1]], he[[1]]$ontology, summary) # all p and FWER = 1 in mole and biol roots
# when skiping works, test if domain=mol works
he2 = go_enrich(ungenes, n_randset=50, domains="molecular_function") ## Error aber warning mit "No GO-annotation for molecular function"; reicht


### erroneous input
# same gene as candidate and background
missass = genes 
missass[1:3,1] = bg_ids[1:3]
res = go_enrich(missass)
# only background defined
res = go_enrich(genes[genes[,2]==0,])
# NA in input
na1 = genes
na1[3,1] = NA
go_enrich(na1)
na2 = genes
na2[3,2] = NA
go_enrich(na2)

##### wilcoxon
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
scores = sample(1:30, length(gene_ids))
genes = data.frame(gene_ids, scores)
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)
head(go_willi[[1]])
go_willi[[2]]
### corner cases
# negative values
genes_rev = data.frame(gene_ids, -scores)
go_willi_rev = go_enrich(genes_rev, test='wilcoxon', n_randsets=100)
head(go_willi_rev[[1]])
forward_p_low = go_willi[[1]][match(go_willi_rev[[1]][,2], go_willi[[1]][,2]),'raw_p_low_rank']
all.equal(go_willi_rev[[1]][,'raw_p_high_rank'], forward_p_low)
# only two genes
go_willi = go_enrich(genes[1:2,], test='wilcoxon', n_randsets=100)
head(go_willi[[1]])
go_willi[[2]]
# only one gene
go_willi = go_enrich(genes[3,], test='wilcoxon', n_randsets=100)
# only two scores - ok
scores = sample(1:2, length(gene_ids), replace=TRUE)
genes = data.frame(gene_ids, scores)
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)
head(go_willi[[1]])
# only one score - works, all p and FWER are 1
genes[genes[,2]==1, 2] = 2
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)
head(go_willi[[1]])
# floating point input
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1')
genes = seq(1.1, 1.7, by=0.1)
names(genes) = gene_ids
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)
head(go_willi[[1]])
# multiple assignment of different values
genes = c(1,2,3)
names(genes) = c('NCAPG', 'APOL4', 'APOL4')
go_enrich(genes, test='wilcoxon')
# NA in input

##### binomial
set.seed(123)
high_human_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_human_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
human_counts = c(sample(20:30, length(high_human_genes)), sample(5:15, length(low_human_genes)))
chimp_counts = c(sample(5:15, length(high_human_genes)), sample(20:30, length(low_human_genes)))
genes = data.frame(gene=c(high_human_genes, low_human_genes), human_counts, chimp_counts)
go_binom = go_enrich(genes, test='binomial', n_randsets=100)
### corner cases
# only one genes - works, p ok, all FWER are 1
go_binom_one = go_enrich(genes[1,], test='binomial', n_randsets=10)
# 0 counts
go_binom_z1 = go_enrich(genes=data.frame(a='G6PD',b=0,c=10), test='binomial', n_randsets=10) # geht
go_binom_z2 = go_enrich(genes=data.frame(a='G6PD',b=10,c=0), test='binomial', n_randsets=10) # geht
# multiple assignment of same values
multi_ok = go_enrich(genes=data.frame(a='G6PD',b=0,c=10)[c(1,1),], test='binomial', n_randsets=10)
## erroneous
go_binom_z3 = go_enrich(genes=data.frame(a='G6PD',b=0,c=0), test='binomial', n_randsets=10)
# negative values
neg_genes = data.frame(gene=c(high_human_genes, low_human_genes), -chimp_counts, -human_counts)
go_binom = go_enrich(neg_genes, test='binomial', n_randsets=10)
# multiple assignment of different values
multi = data.frame(a=c('G6PD','G6PD') ,b=c(0,1),d=c(10,8))
go_enrich(multi, test='binomial')


##### contingency
#func_2x2contingency needs four values per gene. The order of the values are divergence_synonymous divergence_nonsynonymous diversity_syn diversity_nonsyn.
require(GOfuncR)
set.seed(123)
high_substi_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_substi_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2', 'C21orf59')
subs_syn = sample(45:55, length(c(high_substi_genes, low_substi_genes)), replace=T)
subs_non_syn = c(sample(15:25, length(high_substi_genes), replace=T), sample(0:10, length(low_substi_genes)))
vari_syn = sample(25:35, length(c(high_substi_genes, low_substi_genes)), replace=T)
vari_non_syn = c(sample(0:10, length(high_substi_genes), replace=T), sample(10:20, length(low_substi_genes)))
genes = data.frame(genes=c(high_substi_genes, low_substi_genes), vari_syn, vari_non_syn, subs_syn, subs_non_syn)
conti_res = go_enrich(genes, test='contingency', n_randset=100)
### corner cases
# only one gene in one root (C21orf59) - skip other root nodes
conti_root = go_enrich(genes[nrow(genes),], test='contingency', n_randset=100)
# only one gene
conti_res1 = go_enrich(genes[3,], test='contingency', n_randset=100)
# two genes
conti_res2 = go_enrich(genes[3:4,], test='contingency', n_randset=100)
# 0 counts
go_conti_z1 = go_enrich(genes=data.frame(a='G6PD',b=0,c=10,d=0,e=0), test='contingency', n_randsets=10) # geht
go_conti_z2 = go_enrich(genes=data.frame(a='G6PD',b=0,c=0,d=0,e=10), test='contingency', n_randsets=10) # geht
## erroneous
go_conti_z3 = go_enrich(genes=data.frame(a='G6PD',b=0,c=0,d=0,e=0), test='contingency', n_randsets=10) # geht

# NA in input
genesna = genes
genesna[5,4] = NA
go_enrich(genesna, test='contingency')





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
background = read.table('~/R_packages/GOfuncR_package/Ben_background_regions.bed')
candidate = read.table('~/R_packages/GOfuncR_package/Ben_neandertal_deserts.bed')
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




