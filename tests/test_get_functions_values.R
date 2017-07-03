
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
	gos = c('GO:0072025','GO:0072221','GO:0072205','GO:0072235')
	genes = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
	 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
	test_results[["anno_indi"]] = get_anno_genes(go_ids=gos, genes=genes)
	test_results[["anno_indi_all"]] = get_anno_genes(go_ids=gos)
	test_results[["anno_indi_hg20"]] = get_anno_genes(go_ids=gos, genes=genes, ref_genome="grch38")
	test_results[["anno_indi_mus"]] = get_anno_genes(go_ids=gos, ref_genome="grcm38")
	test_results[["anno_indi_mus"]] = get_anno_genes(go_ids=c('GO:007','GO:0072221'), ref_genome="grcm38")

	##### get_anno_categories	
	test_results[["anno_cat"]] = get_anno_categories(c('BTC', 'SPAG5'))
	test_results[["anno_cat_one_wrong"]] = get_anno_categories(c('BTC', '123') , ref_genome='grch38')
	test_results[["anno_cat_none"]] = get_anno_categories('123' , ref_genome='grch38')
	test_results[["anno_cat_mus"]] = get_anno_categories(c('Mus81', 'Papola'), ref_genome='grcm38')

	return(test_results)
}	
	

