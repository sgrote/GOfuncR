

# return a dataframe with genes annotated to GOs 

# input: genes, optional: ref_genome

# output: gene, go_id, go_name


get_anno_categories = function(genes, ref_genome="grch37"){
	
	## Check input	
	if(is.null(genes)){
		stop("Please provide gene-symbols as input.")
	}
	if (!(ref_genome %in% c("grch37","grch38","grcm38"))){
		stop("Please use 'ref_genome=grch37', 'ref_genome=grch38' or 'ref_genome=grcm38'")
	}

	## find annotated GO-categories
	message(paste("find associated categories using ref_genome ",ref_genome,"...",sep=""))
	# remove obsolete terms
	term = term[term[,5]==0,]
	# load GO-annotation
	go_anno = get(paste("go_anno_", ref_genome, sep=""))
	# restrict to genes in annotations
	go_anno = go_anno[go_anno[,2] %in% genes,]
	# restrict annotations to GOs present in term
	out = go_anno[go_anno[,1] %in% term[,4],]

	if(nrow(out)==0){
		warning("No GO-annotations found for the input genes.")
		return(NULL)
	}

	

	# create output
	colnames(out) = c("go_id", "gene")
	out = out[,c(2,1)]
	# add category name
	out = cbind(out, get_names(out$go_id)[,c("name","root_node")])
 	# sort
 	out = out[order(out$gene, out$root_node, out$go_id),]
 	rownames(out) = 1:nrow(out)


	return(out)
}	
