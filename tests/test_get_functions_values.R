
# sample results for 'get_names', 'get_child_nodes', 'get_parent_nodes', 'get_anno_genes'

get_values = function(){
	
	# store all results in a list
	test_results = list()

	##### get names
	test_results[["get_names"]] = get_names(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
	test_results[["get_names_one"]] = get_names('GO:0051082')

	##### get_child_nodes
	test_results[["get_children"]] = get_child_nodes(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
	test_results[["get_children_one"]] = get_child_nodes('GO:0051082')

	##### get_parent_nodes
	test_results[["get_parents"]] = get_parent_nodes(c('GO:0051082', 'GO:123', 'GO:0042254', 'GO:0000109'))
	test_results[["get_parents_root"]] = get_parent_nodes('GO:0003674') #root


	##### get_anno_genes

	### A) given a go_enrich result for GO-IDs with FWER<0.05 (FWER for overrepresentation / high_rank)

	# standard analysis with 100 randomsets 
	set.seed(123)
	gene_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
	 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
	genes = rep(1, length(gene_ids))
	names(genes) = gene_ids
	go_res = go_enrich(genes, n_randset=100)
	test_results[["anno1"]] = get_anno_genes(go_res)
	test_results[["anno_bg1"]] = get_anno_genes(go_res, background=TRUE)
	test_results[["anno_all1"]] = get_anno_genes(go_res, fwer_threshold=5) # its possible to get all annotations

	# background defined
	set.seed(123)
	candi_gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
	bg_gene_ids = c('FGR', 'NPHP1', 'DRD2', 'ABCC10', 'PTBP2', 'JPH4', 'SMARCC2', 'FN1', 'NODAL', 'CYP1A2', 'ACSS1', 'CDHR1', 'SLC25A36', 'LEPR', 'PRPS2', 'TNFAIP3', 'NKX3-1', 'LPAR2', 'PGAM2', 'GAPDHS')
	genes = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
	names(genes) = c(candi_gene_ids, bg_gene_ids)
	go_res = go_enrich(genes, n_randset=100)
	test_results[["anno2"]] = get_anno_genes(go_res, fwer_threshold=0.7)
	test_results[["anno_bg2"]] = get_anno_genes(go_res, background=TRUE, fwer_threshold=0.7)
	test_results[["anno_all2"]] = get_anno_genes(go_res, background=TRUE, fwer_threshold=5)

	# scores for genes and wilcoxon rank sum test
	set.seed(123)
	gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
	 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
	genes = sample(1:30, length(gene_ids))
	names(genes) = gene_ids
	go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)
	test_results[["anno_willi"]] = get_anno_genes(go_willi, fwer_threshold=0.93)
	
	# regions as input
	set.seed(123)
	genes = c(1, rep(0,6))
	names(genes) = c('8:81000000-83000000', '7:1300000-56800000', '7:74900000-148700000',
	 '8:7400000-44300000', '8:47600000-146300000', '9:0-39200000', '9:69700000-140200000')
	go_region = go_enrich(genes, n_randsets=100)
	go_region_hg20 = go_enrich(genes, n_randsets=100, ref_genome="grch38")
	go_region_mus = go_enrich(genes, n_randsets=100, ref_genome="grcm38")
	test_results[["anno_region"]] = get_anno_genes(go_region, fwer_threshold=0.1)
	test_results[["anno_region_bg"]] = get_anno_genes(go_region, fwer_threshold=0.1, background=TRUE)
	test_results[["anno_region_hg20"]] = get_anno_genes(go_region_hg20, fwer_threshold=0.1)
	test_results[["anno_region_mus"]] = get_anno_genes(go_region_mus, fwer_threshold=0.5)


	### B) given GOs and 'optionally' genes directly
	gos = c('GO:0072025','GO:0072221','GO:0072205','GO:0072235')
	genes = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
	 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
	test_results[["anno_indi"]] = get_anno_genes(go_ids=gos, genes=genes, ref_genome="grch37")
	test_results[["anno_indi_all"]] = get_anno_genes(go_ids=gos, ref_genome="grch37")
	test_results[["anno_indi_hg20"]] = get_anno_genes(go_ids=gos, genes=genes, ref_genome="grch38")
	test_results[["anno_indi_mus"]] = get_anno_genes(go_ids=gos, ref_genome="grcm38")

	return(test_results)
}	
	

