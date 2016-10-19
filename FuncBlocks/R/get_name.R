
## return the names of GOs, toghether with root node in a dataframe

get_name=function(go_ids){
	# remove obsolete terms
	term = term[term[,5]==0,]
	# find names for GOs
	out = data.frame(go_ids, term[match(go_ids, term[,4]) ,2:3])
	names(out) = c("go_id", "name", "root_node")
	return(out)
}
