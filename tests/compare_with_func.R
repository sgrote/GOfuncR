
# for each
# * hyper
# * wilcox
# * binomial
# * contingency:
	# run FuncBlocks
	# generate input for FUNC
	# run Func
	# compare output
	# run test for one node in R
	
	
require(FuncBlocks)

# create temp-directory
directory = tempdir()

# plot Func vs FuncBlocks
compi = function(col1, col2, main, xlab, ylab){
	plot(col1, col2, main=main, xlab=xlab, ylab=ylab)
	abline(0,1,col="red")
	comp = all.equal(col1, col2, tolerance=1.5e-6)
	legend("topleft", legend=comp, title="equal", bty="n")
}


### HYPERGEOMETRIC

# run FuncBlocks
candi_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4')
bg_ids = c('C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
scores = c(rep(1,length(candi_ids)), rep(0,length(bg_ids)))
genes_hyper = data.frame(g=c(candi_ids, bg_ids), scores)
set.seed(123)
res_hyper = go_enrich(genes_hyper)
head(res_hyper[[1]])

# prepare input for Func
input_hyper = get_anno_categories(genes_hyper[,1])[,1:2]
input_hyper$score = genes_hyper[match(input_hyper[,1], genes_hyper[,1]),2]
write.table(input_hyper, paste0(directory,"/input_hyper.tsv"), sep="\t", quote=F, row.names=F, col.names=F)

# run Func
system(paste0("func_hyper -i ",directory,"/input_hyper.tsv -t ~/R_packages/FuncBlocks_package/term_tables -o ",directory))
func_weird_out = read.table(paste0(directory,"/groups.txt"), header=T, comment.char="", sep="\t", row.names=NULL)
func_out_hyper = func_weird_out[,1:13]
colnames(func_out_hyper) = colnames(func_weird_out)[2:14]

# compare results
func_out_hyper = func_out_hyper[order(func_out_hyper$FWER_overrepresentation),]
head(func_out_hyper)
head(res_hyper[[1]])
merge_out = merge(func_out_hyper[,c(3,7,8:11)], res_hyper[[1]][,c(2:7,9)])
par(mfrow=c(2,2), pty="s")
compi(merge_out[,3], merge_out[,8], "raw_p_underrep", "Func", "FuncBlocks")
compi(merge_out[,4], merge_out[,9], "raw_p_overrep", "Func", "FuncBlocks")
compi(merge_out[,5], merge_out[,10], "FWER_underrep", "Func", "FuncBlocks")
compi(merge_out[,6], merge_out[,11], "FWER_overrep", "Func", "FuncBlocks")

# run test for one node in R
stich = func_out_hyper[func_out_hyper$node_id=="GO:0044422",]
candi_node = stich$X.genes_with_variable.1_in_node
candi_root = stich$X.genes_with_variable.1_in_root_node
bg_root = stich$X.genes_in_root_node - candi_root
total_node = stich$X.genes_in_node
# under-representation
stich$raw_p_underrepresentation_of_variable.1
phyper(candi_node, candi_root, bg_root, total_node)
# over-respresentation
stich$raw_p_overrepresentation_of_variable.1
phyper(candi_node-1, candi_root, bg_root, total_node, lower.tail=F)



### WILCOXON

# run FuncBlocks
set.seed(123)
high_score_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_score_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
scores = c(sample(20:30, length(high_score_genes)), sample(5:15, length(low_score_genes)))
genes_wilcox = data.frame(genes_wilcox=c(high_score_genes, low_score_genes), scores)
res_wilcox = go_enrich(genes_wilcox, test='wilcoxon')

# prepare input for Func
input_wilcox = get_anno_categories(genes_wilcox[,1])[,1:2]
input_wilcox$score = genes_wilcox[match(input_wilcox[,1], genes_wilcox[,1]),2]
write.table(input_wilcox, paste0(directory,"/input_wilcox.tsv"), sep="\t", quote=F, row.names=F, col.names=F)

# run Func
system(paste0("func_wilcoxon -i ",directory,"/input_wilcox.tsv -t ~/R_packages/FuncBlocks_package/term_tables -o ",directory))
func_weird_out = read.table(paste0(directory,"/groups.txt"), header=T, comment.char="", sep="\t", row.names=NULL)
func_out_wilcox = func_weird_out[,1:12]
colnames(func_out_wilcox) = colnames(func_weird_out)[2:13]

# compare results
func_out_wilcox = func_out_wilcox[order(func_out_wilcox$FWER_high),]
head(func_out_wilcox)
head(res_wilcox[[1]])
merge_out = merge(func_out_wilcox[,c(3,7:10)], res_wilcox[[1]][,c(2:7,9)])
par(mfrow=c(2,2), pty="s")
compi(merge_out[,2], merge_out[,7], "raw_p_low_rank", "Func", "FuncBlocks")
compi(merge_out[,3], merge_out[,8], "raw_p_high_rank", "Func", "FuncBlocks")
compi(merge_out[,4], merge_out[,9], "FWER_low_rank", "Func", "FuncBlocks")
compi(merge_out[,5], merge_out[,10], "FWER_high_rank", "Func", "FuncBlocks")

# run test for one node in R

# a) as implemented in Func and FuncBlocks
stich = func_out_wilcox[func_out_wilcox$node_id=="GO:0008152",]
ranksum = stich$sum_of_ranks_in_node
sum_nties = 0
n = stich$X.genes_in_node
N = stich$X.genes_outside_node
C = ranksum - ((n*(n+1))/2)
z = C - (n*N*0.5)
sigma = sqrt((n*N/12) * ((n+N+1)-sum_nties / ((n+N)*(N+n-1.))))
# p low rank
corr = -0.5
stich$raw_p_low_ranks
pnorm( (z-corr) / sigma, 0, 1, 1, 0 )
# p high rank
corr = 0.5
stich$raw_p_high_ranks
1 - pnorm( (z-corr) / sigma, 0, 1, 1, 0 )

# b) ordinary wilcox test - not exacactly the same (different correction?)
anno_node = get_anno_genes("GO:0008152", genes=genes_wilcox[,1])
anno_root = get_anno_genes(get_ids(get_names("GO:0008152")[,3])[,3], genes=genes_wilcox[,1]) 
anno_root = anno_root[!(anno_root[,2] %in% anno_node[,2]),] # genes outside node; not genes in root
anno_node$score = genes_wilcox[match(anno_node[,2], genes_wilcox[,1]),2]
anno_root$score = genes_wilcox[match(anno_root[,2], genes_wilcox[,1]),2]
# p low rank
stich$raw_p_low_ranks
wilcox.test(anno_node$score, anno_root$score, alternative="less")
# p high rank
stich$raw_p_high_ranks
wilcox.test(anno_node$score, anno_root$score, alternative="greater")



### BINOMIAL

# run FuncBlocks
set.seed(123)
high_human_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_human_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
human_counts = c(sample(20:30, length(high_human_genes)), sample(5:15, length(low_human_genes)))
chimp_counts = c(sample(5:15, length(high_human_genes)), sample(20:30, length(low_human_genes)))
genes_binom = data.frame(gene=c(high_human_genes, low_human_genes), chimp_counts, human_counts)
res_binom = go_enrich(genes_binom, test='binomial', silent=T)

# prepare input for Func
input_binom = get_anno_categories(genes_binom[,1])[,1:2]
input_binom = cbind(input_binom, genes_binom[match(input_binom[,1], genes_binom[,1]),2:3])
write.table(input_binom, paste0(directory,"/input_binom.tsv"), sep="\t", quote=F, row.names=F, col.names=F)

# run Func
system(paste0("func_binom -i ",directory,"/input_binom.tsv -t ~/R_packages/FuncBlocks_package/term_tables -o ",directory))
func_weird_out = read.table(paste0(directory,"/groups.txt"), header=T, comment.char="", sep="\t", row.names=NULL)
func_out_binom = func_weird_out[,1:11]
colnames(func_out_binom) = colnames(func_weird_out)[2:12]

# compare results
func_out_binom = func_out_binom[order(func_out_binom$FWER_2),]
head(func_out_binom)
head(res_binom[[1]])
merge_out = merge(func_out_binom[,c(3,6:9)], res_binom[[1]][,c(2:7)])
par(mfrow=c(2,2), pty="s")
compi(merge_out[,2], merge_out[,7], "raw_p_high_A", "Func", "FuncBlocks")
compi(merge_out[,3], merge_out[,8], "raw_p_high_B", "Func", "FuncBlocks")
compi(merge_out[,4], merge_out[,9], "FWER_high_A", "Func", "FuncBlocks")
compi(merge_out[,5], merge_out[,10], "FWER_high_B", "Func", "FuncBlocks")

# run test for one node in R
stich = func_out_binom[func_out_binom$node_id=="GO:0044422",]
root = func_out_binom[func_out_binom$node_name=="cellular_component",]
a_node = stich$sum_gene_associated_variable_1
b_node = stich$sum_gene_associated_variable_2
a_root = root$sum_gene_associated_variable_1
b_root = root$sum_gene_associated_variable_2
p_a = a_root / (a_root + b_root)
p_b = b_root / (a_root + b_root)
# high_A
stich$raw_p_overrepresentation_1
binom.test(c(a_node, b_node), p = p_a, alternative="greater")
binom.test(c(b_node, a_node), p = p_b, alternative="less") # same
# high_B
stich$raw_p_overrepresentation_2
binom.test(c(b_node, a_node), p = p_b, alternative="greater")



### CONTINGENCY

# run FuncBlocks
set.seed(123)
high_substi_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_substi_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2', 'C21orf59')
subs_syn = sample(45:55, length(c(high_substi_genes, low_substi_genes)), replace=T)
subs_non_syn = c(sample(15:25, length(high_substi_genes), replace=T), sample(0:10, length(low_substi_genes)))
vari_syn = sample(25:35, length(c(high_substi_genes, low_substi_genes)), replace=T)
vari_non_syn = c(sample(0:10, length(high_substi_genes), replace=T), sample(10:20, length(low_substi_genes)))
genes_conti = data.frame(genes=c(high_substi_genes, low_substi_genes), vari_syn, vari_non_syn, subs_syn, subs_non_syn)
res_conti = go_enrich(genes_conti, test='contingency', silent=T)

# prepare input for Func
input_conti = get_anno_categories(genes_conti[,1])[,1:2]
input_conti = cbind(input_conti, genes_conti[match(input_conti[,1], genes_conti[,1]),2:5])
write.table(input_conti, paste0(directory,"/input_conti.tsv"), sep="\t", quote=F, row.names=F, col.names=F)

# run Func
system(paste0("func_2x2contingency -i ",directory,"/input_conti.tsv -t ~/R_packages/FuncBlocks_package/term_tables -o ",directory))
func_weird_out = read.table(paste0(directory,"/groups.txt"), header=T, comment.char="", sep="\t", row.names=NULL)
func_out_conti = func_weird_out[,1:13]
colnames(func_out_conti) = colnames(func_weird_out)[2:14]

# compare results ## TODO: hier auf eine neue Version von Func warten
# (aber vor bug-fix in FuncBlocks haben die Ergebnisse auch uebereingestimmt)
func_out_conti = func_out_conti[order(func_out_conti$FWER_3.4),]
head(func_out_conti)
head(res_conti[[1]])
merge_out = merge(func_out_conti[,c(3,8:11)], res_conti[[1]][,c(2:7)])
par(mfrow=c(2,2), pty="s")
compi(merge_out[,2], merge_out[,7], "raw_p_high_C_D", "Func", "FuncBlocks")
compi(merge_out[,3], merge_out[,8], "raw_p_high_A_B", "Func", "FuncBlocks")
compi(merge_out[,4], merge_out[,9], "FWER_high_C_D", "Func", "FuncBlocks")
compi(merge_out[,5], merge_out[,10], "FWER_high_C_D", "Func", "FuncBlocks")

# run test for one node in R
# (two stichs because one of the p-vals is always 1) 
conti = res_conti[[1]]
# apparently a chisq-test is performed, and if A/B > C/D chisq p-val is assigned high_A/B, else high_C/D and the other one is set to 1
# p_high_AD
stich1 = conti[conti$node_id=="GO:0044422",]
stich1_anno = get_anno_genes(stich1$node_id, genes=genes_conti[,1])
stich1_anno = cbind(stich1_anno, genes_conti[match(stich1_anno[,2], genes_conti[,1]),2:5])
conti_tab = matrix(colSums(stich1_anno[,3:6]), ncol=2) # matrix is filled col-wise
chisq.test(conti_tab, correct=FALSE)
stich1$raw_p_high_AB
# p_high_CD
stich2 = conti[conti$node_id=="GO:0046903",]
stich2_anno = get_anno_genes(stich2$node_id, genes=genes_conti[,1])
stich2_anno = cbind(stich2_anno, genes_conti[match(stich2_anno[,2], genes_conti[,1]),2:5])
conti_tab2 = matrix(colSums(stich2_anno[,3:6]), ncol=2) # matrix is filled col-wise
chisq.test(conti_tab2, correct=FALSE)
stich2$raw_p_high_CD








