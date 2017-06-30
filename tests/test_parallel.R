library("FuncBlocks")
library("parallel")

message("test parallel...")

candi1_gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'NPHP1', 'DRD2', 'FN1', 'NODAL')
candi2_gene_ids = c('C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'ABCC10', 'PTBP2', 'CYP1A2', 'ACSS1')
candi3_gene_ids = c('BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2', 'JPH4', 'SMARCC2', 'CDHR1', 'SLC25A36')
bg_gene_ids = c('FGR', 'LEPR', 'PRPS2', 'TNFAIP3', 'NKX3-1', 'LPAR2', 'PGAM2', 'GAPDHS')
input=list()


##### a) hypergeometric

input[["genes1"]] = c(rep(1,length(candi1_gene_ids)), rep(0,length(bg_gene_ids)))
names(input[["genes1"]]) = c(candi1_gene_ids, bg_gene_ids)
input[["genes2"]] = c(rep(1,length(candi2_gene_ids)), rep(0,length(bg_gene_ids)))
names(input[["genes2"]]) = c(candi2_gene_ids, bg_gene_ids)
input[["genes3"]] = c(rep(1,length(candi3_gene_ids)), rep(0,length(bg_gene_ids)))
names(input[["genes3"]]) = c(candi3_gene_ids, bg_gene_ids)

# a) loop
res = list()
for(i in 1:3){
	set.seed(123)
	res[[i]] = go_enrich(input[[i]], n_randset=50, silent=T)
}

# b) parallel 
parares = mclapply(1:3, function(x){
	set.seed(123)
	go_enrich(input[[x]], n_randset=50, silent=T)
})

# all 3 parares results are like res[[3]] # does not work with FuncBlocks_1.2.3
for (i in 1:3){
	if(!(isTRUE(all.equal(res[[i]], parares[[i]])))){
		stop("Error in test_parallel.R (hypergeometric)", silent=T)
	}
}


##### b) wilcoxon

input[["genes1"]] = sample(1:30, length(candi1_gene_ids) + length(bg_gene_ids))
names(input[["genes1"]]) = c(candi1_gene_ids, bg_gene_ids)
input[["genes2"]] = sample(1:30, length(candi2_gene_ids) + length(bg_gene_ids))
names(input[["genes2"]]) = c(candi2_gene_ids, bg_gene_ids)
input[["genes3"]] = sample(1:30, length(candi3_gene_ids) + length(bg_gene_ids))
names(input[["genes3"]]) = c(candi3_gene_ids, bg_gene_ids)

# a) loop
res = list()
for(i in 1:3){
	set.seed(123)
	res[[i]] = go_enrich(input[[i]], n_randset=50, test="wilcoxon", silent=T)
}

# b) parallel
parares = mclapply(1:3, function(x){
	set.seed(123)
	go_enrich(input[[x]], n_randset=50, test="wilcoxon", silent=T)
})

# all 3 parares results are like res[[3]] # does not work with FuncBlocks_1.2.3
for (i in 1:3){
	if(!(isTRUE(all.equal(res[[i]], parares[[i]])))){
		stop("Error in test_parallel.R (wilcoxon)")
	}
}
