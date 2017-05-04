
go_enrich_values = function(){
	set.seed(123)

	# store all results in a list
	test_results = list()

	###########################################################################################

	##### standard parameters
	set.seed(123)
	gene_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
#	gene_ids = c('NCAPG', 'QUATSCH1', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
	genes = rep(1, length(gene_ids))
	names(genes) = gene_ids
	test_results[["res1"]] = go_enrich(genes, n_randset=100)
	test_results[["res1_hg20"]] = go_enrich(genes, n_randset=100, ref_genome='grch38')
	### corner cases
	# one gene
	test_results[["res2"]] = go_enrich(genes[1], n_randsets=100)
	### mouse
	gene_ids = c('Arsi', 'Mapk4', 'Papola', 'Tfrc', 'Bak1', 'Fopnl', 'Mus81', 'Opa3', 'Npcd')
	genes = rep(1, length(gene_ids))
	names(genes) = gene_ids
	test_results[["res_mouse"]] = go_enrich(genes, n_randsets=100, ref_genome='grcm38')


	##### standard parameters - background defined
	set.seed(123)
	candi_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4')
	bg_ids = c('C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
	genes = c(rep(1,length(candi_ids)), rep(0,length(bg_ids)))
	names(genes) = c(candi_ids, bg_ids)
	test_results[["res3"]] = go_enrich(genes, n_randsets=100)
	### corner cases
	# one background gene
	test_results[["res4"]] = go_enrich(genes[1:(length(candi_ids)+1)], n_randsets=100)
	# one candidate gene
	test_results[["res5"]] = go_enrich(genes[c(1,(length(candi_ids)+1):length(genes))], n_randsets=100)


	##### wilcoxon
	set.seed(123)
	gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
	genes = sample(1:30, length(gene_ids))
	names(genes) = gene_ids
	test_results[["go_willi1"]] = go_enrich(genes, test='wilcoxon', n_randsets=100)
	### corner cases
	# negative values
	genes_rev=-genes
	test_results[["go_willi2"]] = go_enrich(genes_rev, test='wilcoxon', n_randsets=100)
	# only two genes
	test_results[["go_willi3"]] = go_enrich(genes[1:2], test='wilcoxon', n_randsets=100)
	# only two scores
	genes = sample(1:2, length(gene_ids), replace=TRUE)
	names(genes) = gene_ids
	test_results[["go_willi4"]] = go_enrich(genes, test='wilcoxon', n_randsets=100)
	# only one score - works, all p and FWER are 1
	genes[genes==1] = 2
	test_results[["go_willi5"]] = go_enrich(genes, test='wilcoxon', n_randsets=10)
	# floating point input
	gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1')
	genes = seq(1.1, 1.7, by=0.1)
	names(genes) = gene_ids
	test_results[["go_willi6"]] = go_enrich(genes, test='wilcoxon', n_randsets=100)


	##### gene_len
	set.seed(123)
	gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2')
	genes = rep(1, length(gene_ids))
	names(genes) = gene_ids
	test_results[["res_len"]] = go_enrich(genes, gene_len=TRUE, n_randset=100)
	test_results[["res_len_hg20"]] = go_enrich(genes, gene_len=TRUE, n_randset=100, ref_genome='grch38')
	### mouse
	gene_ids = c('Arsi', 'Mapk4', 'Papola', 'Tfrc', 'Bak1', 'Fopnl', 'Mus81', 'Opa3', 'Npcd')
	genes = rep(1, length(gene_ids))
	names(genes) = gene_ids
	test_results[["res_mouse_len"]] = go_enrich(genes, gene_len=TRUE, n_randsets=100, ref_genome='grcm38')


	##### genomic regions
	set.seed(123)
	genes = c(1,1, rep(0,6))
	names(genes) = c('8:81000000-83000000', '3:76500000-90500000', '7:1300000-56800000', '7:74900000-148700000',
	 '8:7400000-44300000', '8:47600000-146300000', '9:0-39200000', '9:69700000-140200000')
	test_results[["go_region"]] = go_enrich(genes, n_randsets=100)
	### hg20
	test_results[["go_region_hg20"]] = go_enrich(genes, n_randsets=100, ref_genome='grch38')
	### mouse
	test_results[["go_region_mouse"]] = go_enrich(genes, n_randsets=100, ref_genome='grcm38')

	### circ_chrom
	set.seed(123)
	test_results[["go_circ1"]] = go_enrich(genes[c(1,5,6)], n_randsets=100, circ_chrom=TRUE)
	test_results[["go_circ1_hg20"]] = go_enrich(genes[c(1,5,6)], n_randsets=100, circ_chrom=TRUE, ref_genome='grch38')
	# warning about unused chromosomes
	test_results[["go_circ2"]] = go_enrich(genes[-2], n_randsets=100, circ_chrom=TRUE)

	return(test_results)
}


