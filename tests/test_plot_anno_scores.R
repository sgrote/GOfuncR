
# run in GOfuncR/tests

library(GOfuncR)

pdf("test_plot_anno_scores.pdf", width=10)

	################## (1) hyper

	set.seed(123)

	## run Func
	gene_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
	 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
	genes = rep(1, length(gene_ids))
	names(genes) = gene_ids
	go_res = go_enrich(genes, n_randsets=20, silent=T)

	### run Func - defined background
	candi_gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
	bg_gene_ids = c('FGR', 'NPHP1', 'DRD2', 'ABCC10', 'PTBP2', 'JPH4', 'SMARCC2', 'FN1', 'NODAL', 'CYP1A2', 'ACSS1', 'CDHR1', 'SLC25A36', 'LEPR', 'PRPS2', 'TNFAIP3', 'NKX3-1', 'LPAR2', 'PGAM2', 'GAPDHS')
	genes = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
	names(genes) = c(candi_gene_ids, bg_gene_ids)
	go_bg = go_enrich(genes, n_randsets=20, silent=T)
	### skip one root
	set.seed(123)
	go_bg_skip = go_enrich(genes, n_randsets=20, domains=c('molecular_function', 'cellular_component'), silent=T)
	go_bg_skip2 = go_enrich(genes, n_randsets=20, domains=c('molecular_function', 'biological_process'), silent=T)
	go_bg_skip3 = go_enrich(genes, n_randsets=20, domains=c('cellular_component'), silent=T)


	### plot 
	# normal input (check order-preserve)
	x = plot_anno_scores(go_res, go_ids=c('GO:0072025','GO:0072221','GO:0072235', 'GO:0044765'))
	print(x)
	plot_anno_scores(go_bg, go_ids=c('GO:0005634','GO:0004945','GO:0008289','GO:0005737','GO:0071495'))
	# corner case
	plot_anno_scores(go_bg, go_ids=c('GO:0005623')) # only one GO
	plot_anno_scores(go_bg, go_ids=c('GO:0000166', 'GO:0000287', 'GO:0000981')) # no candidate annotated
	# skipped root node
	plot_anno_scores(go_bg_skip, go_bg_skip[[1]][1:5,2]) # no biol_process
	plot_anno_scores(go_bg_skip2, go_bg_skip2[[1]][1:5,2]) # no cell_comp
	plot_anno_scores(go_bg_skip3, go_bg_skip3[[1]][1:5,2]) # only cell_comp



	################## (2) wilcoxon
	library(GOfuncR)
	set.seed(123)
	high_score_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
	low_score_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
	genes = c(sample(20:30, length(high_score_genes)), sample(5:15, length(low_score_genes)))
	names(genes) = c(high_score_genes, low_score_genes)
	go_willi = go_enrich(genes, test='wilcoxon', n_randsets=20, silent=T)
	go_willi_skip = go_enrich(genes, test='wilcoxon', n_randsets=20, domains=c('cellular_component'), silent=T)

	## plot
	# normal input
	w = plot_anno_scores(go_willi, go_ids=c('GO:0072025','GO:0072221','GO:0072235', 'GO:0044765'))
	print(w)
	plot_anno_scores(go_willi, go_ids=c('GO:0005634','GO:0004945','GO:0008289','GO:0005737','GO:0071495'))
	# corner cases
	plot_anno_scores(go_willi, go_ids=c('GO:0005634')) # only one GO
	plot_anno_scores(go_willi_skip, go_ids=go_willi_skip[[1]][1:5,2]) # skipped root node

	# mouse
	set.seed(123)
	high_score_genes = c('G6pd', 'Gck', 'Gys1', 'Hk2', 'Pygl', 'Slc2a8', 'Ugp2', 'Zwint', 'Engase')
	low_score_genes = c('Cacng2', 'Agtr1', 'Ano1', 'Btbd3', 'Mtus1', 'Calb1', 'Gyg1', 'Pax2')
	genes = c(sample(20:30, length(high_score_genes)), sample(5:15, length(low_score_genes)))
	names(genes) = c(high_score_genes, low_score_genes)
	go_willi_mus = go_enrich(genes, test='wilcoxon', n_randsets=20, ref_genome='grcm38', silent=T)
	plot_anno_scores(go_willi_mus, go_ids=c('GO:0072025','GO:0072221','GO:0072235', 'GO:0044765')) # check mouse genome is used in get_anno_genes (prints to console 'using ref-genome...')


	####################  (3) binomial
	require(GOfuncR)
	set.seed(123)
	high_human_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
	low_human_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
	human_counts = c(sample(20:30, length(high_human_genes)), sample(5:15, length(low_human_genes)))
	chimp_counts = c(sample(5:15, length(high_human_genes)), sample(20:30, length(low_human_genes)))
	genes = data.frame(gene=c(high_human_genes, low_human_genes), chimp_counts, human_counts)
	go_binom = go_enrich(genes, test='binomial', n_randsets=20, silent=T)
	go_binom_skip = go_enrich(genes, test='binomial', n_randsets=20, domains=c('molecular_function', 'cellular_component'), silent=T)
	go_binom_skip2 = go_enrich(genes, test='binomial', n_randsets=20, domains=c('molecular_function', 'biological_process'), silent=T)
	go_binom_skip3 = go_enrich(genes, test='binomial', n_randsets=20, domains=c('cellular_component'), silent=T)
	### new
	b = plot_anno_scores(go_binom, go_ids = head(go_binom[[1]])$node_id)
	print(b)
	plot_anno_scores(go_binom, go_ids = tail(go_binom[[1]])$node_id)
	plot_anno_scores(go_binom_skip, go_ids = head(go_binom_skip[[1]])$node_id[1:5])
	plot_anno_scores(go_binom_skip2, go_ids = head(go_binom_skip2[[1]])$node_id[1:3])
	# only one node
	plot_anno_scores(go_binom_skip3, go_ids = head(go_binom_skip3[[1]])$node_id[4])


	#################### (4) contingency
	# 2x2 contingency
	require(GOfuncR)
	set.seed(123)
	high_substi_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
	low_substi_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2', 'C21orf59')
	subs_syn = sample(45:55, length(c(high_substi_genes, low_substi_genes)), replace=T)
	subs_non_syn = c(sample(15:25, length(high_substi_genes), replace=T), sample(0:10, length(low_substi_genes)))
	vari_syn = sample(25:35, length(c(high_substi_genes, low_substi_genes)), replace=T)
	vari_non_syn = c(sample(0:10, length(high_substi_genes), replace=T), sample(10:20, length(low_substi_genes)))
	genes = data.frame(genes=c(high_substi_genes, low_substi_genes), vari_syn, vari_non_syn, subs_syn, subs_non_syn)
	go_conti = go_enrich(genes, test="contingency", n_randsets=20, silent=T)

	## plot
	x = plot_anno_scores(go_conti, go_ids=c('GO:0072025','GO:0072221','GO:0072235', 'GO:0044765'))
	print(x)
	# only one node
	plot_anno_scores(go_conti, go_ids=c('GO:0072025'))


dev.off()

#### errors
tryCatch(plot_anno_scores("bla"), error=function(cond) print(cond))

# no genes annotated 
# (plot_anno_scores checks if go_ids in res[[1]], nodes without annotations are included in res in the first place)
tryCatch(plot_anno_scores(go_bg, go_ids=c('GO:0000009', 'GO:0000010', 'GO:0000014')), error=function(cond) print(cond))
# GO not in res-root-nodes (also skipped domain, should be the same error as above)
tryCatch(plot_anno_scores(go_bg_skip3, go_ids=c('GO:0005634','GO:0004945')), error=function(cond) print(cond))

message(paste0(getwd(),"/test_plot_anno_scores.pdf created."))
