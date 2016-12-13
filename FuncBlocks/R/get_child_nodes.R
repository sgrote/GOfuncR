# return a dataframe with parent and child nodes and distances 

get_child_nodes = function(go_ids){

	go_ids = as.character(go_ids)
	# remove obsolete terms
	term = term[term[,5]==0,]
	# get term-IDs of GOs
	go_ids_term = term[term[,4] %in% go_ids, 1]
	if(length(go_ids_term) == 0){
		stop("GO-IDs not found in ontology or obsolete.")
	}
	# go-categories that are children of *go* (cols graph_path: id, parent, child, relation, distance, ...)
	# term-IDs and distance
	go_children_id = graph_path[graph_path[,2] %in% go_ids_term,]
	# names and GO-IDS
	parent = term[match(go_children_id[,2],term[,1]) ,4]
	# GO-IDs of child nodes (and the node itself)
	children = term[match(go_children_id[,3],term[,1]) ,4]
	out = data.frame(parent, children, go_children_id[,5])
	colnames(out) = c("parent_go_id","child_go_id","distance")
	# sort
	out = out[order(out[,1],out[,3],out[,2]),]

	return(out)
}

# TODO: add example with adding name via get_names
