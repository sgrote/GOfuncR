
## return the IDs and names of GO-categories that match ONE name

get_ids = function(go_name){
	if(length(go_name)!=1){
		stop("Please use a single GO-category name.")
	}
	# restrict to 'real' GO-IDs
	term = term[substr(term[,4],1,3)=="GO:",]
	# remove obsolete
	term = term[term[,5]==0,]
	# grep for name
	out = term[grep(go_name, term[,2], ignore.case=TRUE),2:4]
	if(nrow(out)==0){
		stop("No matches found.")
	}
	out = out[order(out[,2], out[,3]),]
	rownames(out) = 1:nrow(out)
	colnames(out) = c("node_name", "root_node", "go_id")
	
	return(out)
}

