
# for each
# * hyper
# * wilcox
# * binomial:
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
genes = data.frame(g=c(candi_ids, bg_ids), scores)
set.seed(123)
res = go_enrich(genes)
head(res[[1]])

# prepare input for Func
input = get_anno_categories(genes[,1])[,1:2]
input$score = genes[match(input[,1], genes[,1]),2]
write.table(input, paste0(directory,"/input_hyper.tsv"), sep="\t", quote=F, row.names=F, col.names=F)

# run Func
system(paste0("func_hyper -i ",directory,"/input_hyper.tsv -t ~/R_packages/FuncBlocks_package/term_tables -o ",directory))
func_weird_out = read.table(paste0(directory,"/groups.txt"), header=T, comment.char="", sep="\t", row.names=NULL)
func_out = func_weird_out[,1:13]
colnames(func_out) = colnames(func_weird_out)[2:14]

# compare results
func_out = func_out[order(func_out$FWER_overrepresentation),]
head(func_out)
head(res[[1]])
merge_out = merge(func_out[,c(3,7,8:11)], res[[1]][,c(2:7,9)])
par(mfrow=c(2,2), pty="s")
compi(merge_out[,3], merge_out[,8], "raw_p_underrep", "Func", "FuncBlocks")
compi(merge_out[,4], merge_out[,9], "raw_p_overrep", "Func", "FuncBlocks")
compi(merge_out[,5], merge_out[,10], "FWER_underrep", "Func", "FuncBlocks")
compi(merge_out[,6], merge_out[,11], "FWER_overrep", "Func", "FuncBlocks")

# run test for one node in R
stich = func_out[func_out$node_id=="GO:0044422",]
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
genes = data.frame(genes=c(high_score_genes, low_score_genes), scores)
res = go_enrich(genes, test='wilcoxon')

# prepare input for Func
input = get_anno_categories(genes[,1])[,1:2]
input$score = genes[match(input[,1], genes[,1]),2]
write.table(input, paste0(directory,"/input_wilcox.tsv"), sep="\t", quote=F, row.names=F, col.names=F)

# run Func
system(paste0("func_wilcoxon -i ",directory,"/input_wilcox.tsv -t ~/R_packages/FuncBlocks_package/term_tables -o ",directory))
func_weird_out = read.table(paste0(directory,"/groups.txt"), header=T, comment.char="", sep="\t", row.names=NULL)
func_out = func_weird_out[,1:12]
colnames(func_out) = colnames(func_weird_out)[2:13]

# compare results
func_out = func_out[order(func_out$FWER_high),]
head(func_out)
head(res[[1]])
merge_out = merge(func_out[,c(3,7:10)], res[[1]][,c(2:7,9)])
par(mfrow=c(2,2), pty="s")
compi(merge_out[,2], merge_out[,7], "raw_p_low_rank", "Func", "FuncBlocks")
compi(merge_out[,3], merge_out[,8], "raw_p_high_rank", "Func", "FuncBlocks")
compi(merge_out[,4], merge_out[,9], "FWER_low_rank", "Func", "FuncBlocks")
compi(merge_out[,5], merge_out[,10], "FWER_high_rank", "Func", "FuncBlocks")

# run test for one node in R

# a) as implemented in Func and FuncBlocks
stich = func_out[func_out$node_id=="GO:0008152",]
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
anno_node = get_anno_genes("GO:0008152", genes=genes[,1])
anno_root = get_anno_genes(get_ids(get_names("GO:0008152")[,3])[,3], genes=genes[,1]) 
anno_root = anno_root[!(anno_root[,2] %in% anno_node[,2]),] # genes outside node; not genes in root
anno_node$score = genes[match(anno_node[,2], genes[,1]),2]
anno_root$score = genes[match(anno_root[,2], genes[,1]),2]
# p low rank
stich$raw_p_low_ranks
wilcox.test(anno_node$score, anno_root$score, alternative="less")
# p high rank
stich$raw_p_high_ranks
wilcox.test(anno_node$score, anno_root$score, alternative="greater")







