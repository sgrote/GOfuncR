# return a dataframe with child and parent nodes and distances 

get_parent_nodes = function(go_ids){

	go_ids = as.character(go_ids)
	# remove obsolete terms
	term = term[term[,5]==0,]
	# get term-IDs of GOs
	go_ids_term = term[term[,4] %in% go_ids, 1]
	if(length(go_ids_term) == 0){
		stop("GO-IDs not found in ontology or obsolete.")
	}
	# go-categories that are parents of *go* (cols graph_path: id, parent, child, relation, distance, ...)
	# term-IDs and distance
	go_parent_id = graph_path[graph_path[,3] %in% go_ids_term,]
	# names and GO-IDS
	parent = term[match(go_parent_id[,2],term[,1]) ,4]
	# GO-IDs of child nodes (and the node itself)
	children = term[match(go_parent_id[,3],term[,1]) ,4]
	out = data.frame(children, parent, go_parent_id[,5])
	colnames(out) = c("child_go_id","parent_go_id","distance")
	out[,1:2] = apply(out[,1:2], 2, as.character)
	# sort
	out = out[order(out[,1],out[,3],out[,2]),]
	# remove 'all'-root node
	out = out[out[,2] != "all",]
	
	rownames(out) = 1:nrow(out)
	
	return(out)
}

# TODO: add example with adding name via get_names
