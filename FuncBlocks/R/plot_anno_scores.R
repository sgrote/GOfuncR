# take a result from go_enrich() and a vector of GO-IDs to plot annotated scores
# plotting type depends on the performed test in go_enrich which is automatically recognized
# (fwer_threshold is not supported anymore, to not discriminate high_A or high_B in binomial, contingency)


plot_anno_scores = function(res, go_ids){
	

	### check input
	# check that res could be go_enrich-output
	if (!(is.list(res) && all(names(res) == c("results","genes","ref_genome")))){
		stop("Please use an object returned from go_enrich as input (list with 3 elements).")
	}
	# infer test
	in_genes = res[[2]]
	if (ncol(in_genes) == 2){
		if(all(in_genes[,2] %in% c(1,0))){
			test = "hyper"
		} else {
			test = "wilcoxon"
		}
	} else if(ncol(in_genes) == 3){
		test = "binomial"
	} else if(ncol(in_genes) == 5){
		test = "contingency"
	} else {
		stop("Identification of test failed.")
	}
	# default root-ids for stable colors (TODO:remove default if onto is input); not needed for conti...
	def_root_ids = c("GO:0003674","GO:0005575","GO:0008150")
	pie_cols = c("#F15A60","#7BC36A","#599BD3","#F9A75B","#9E67AB","#CE7058","#D77FB4")
	root_cols = data.frame(root=def_root_ids, col=pie_cols[1:length(def_root_ids)], stringsAsFactors=FALSE)
	
	# check if background genes are defined (optional for hyper)
	if (test == "hyper" & all(in_genes[,2] == 1)){
		bgdef = FALSE
	} else {
		bgdef = TRUE
	}
	# just in case
	go_ids = as.character(go_ids)
	# keep order of input GO's (which gets messed up in get_anno_genes by *apply)
	ordere = data.frame(go_ids, rank=1:length(go_ids))
	
	# get IDs for root_nodes
	root_names = unique(res[[1]][,1])
	root_names = term[match(root_names, term[,2]) ,]
	root_ids = root_names[,4]	
	
	# get annotation for go_ids and root-nodes
	if (bgdef){
		# background defined: restrict to input genes
		anno = get_anno_genes(go_id=go_ids, genes=in_genes[,1], ref_genome=res[[3]])
		root_anno = get_anno_genes(go_id=root_ids, genes=in_genes[,1], ref_genome=res[[3]])
	} else {
		# background not defined: get all genes annotations
		anno = get_anno_genes(go_id=go_ids, ref_genome=res[[3]])
		root_anno = get_anno_genes(go_id=root_ids, ref_genome=res[[3]])
	}
	if (is.null(anno)) return(invisible(anno)) # no annotations - warning from get_anno_genes
	
	if (test != "wilcoxon"){
		# add scores to nodes and root-nodes
		anno_scores = cbind(anno, in_genes[match(anno[,2], in_genes[,1]), 2:ncol(in_genes)])
		root_anno_scores = cbind(root_anno, in_genes[match(root_anno[,2], in_genes[,1]), 2:ncol(in_genes)])
		if (!bgdef){
			anno_scores[is.na(anno_scores[,3]), 3] = 0 # default 0 for hyper
			root_anno_scores[is.na(root_anno_scores[,3]), 3] = 0
		}
		# aggregate scores in nodes and root-nodes
		if (test == "hyper"){ 
			# counts of 1 and 0 genes in a node
			aggrego = tapply(anno_scores[,3], anno_scores[,1], function(x) c(sum(x), length(x)-sum(x)))
			aggrego = data.frame(go_id = names(aggrego), do.call(rbind, aggrego))
			root_aggrego = tapply(root_anno_scores[,3], root_anno_scores[,1], function(x) c(sum(x), length(x)-sum(x)))
			root_aggrego = data.frame(go_id = names(root_aggrego), do.call(rbind, root_aggrego))
		} else { 
			# sums of scores in a node
			aggrego = aggregate(anno_scores[,3:ncol(anno_scores)], list(go_id=anno_scores[,1]), sum)
			root_aggrego = aggregate(root_anno_scores[,3:ncol(root_anno_scores)], list(go_id=root_anno_scores[,1]), sum)
		}
		# add colors and root_node_name
		root_aggrego$root_name = get_names(root_aggrego[,1])[,2]
		root_aggrego$root_col = root_cols[match(root_aggrego[,1], root_cols[,1]), 2]
		# merge nodes with root node info
		matched_root_name = get_names(aggrego[,1])[,3] # get names
		aggrego$root_id = root_names[match(matched_root_name, root_names[,2]), 4]
		aggrego = cbind(aggrego, root_aggrego[match(aggrego$root_id, root_aggrego[,1]), 2:ncol(root_aggrego)])
	}
	# get original order
	aggrego = aggrego[order(ordere[match(aggrego$go_id, ordere$go_ids), 2]),]
	rownames(aggrego) = 1:nrow(aggrego)
	
	# plot and get stats returned
	if (test == "hyper"){
		stats = plot_hyper(aggrego, root_aggrego)
	} else if (test == "binomial"){
		stats = plot_binomial(aggrego, root_aggrego)
	}
	
	# wilcoxon:
		
		# table root_info (go_id -> node-id)
		# (only test where scores are not collapsed, would be too redundant to have everything in one table) 
	
		# create output with list of 3 tables(node-anno, root-anno, node->root-mapping)
		
		# add colors and root_node_name
	
	
	# remove unused columns for output (maybe add GO-name, useful also because plot just has IDs)
	return(invisible(stats))
}
	
	
	

