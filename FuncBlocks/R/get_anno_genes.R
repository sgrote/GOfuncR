

# return a dataframe with genes annotated to (enriched) GOs 

# input: result list produced by go_enrich, optional: FWER, background
# alternative input: go_ids, optional: genes

# output: go_id, gene, (FWER, candiate/score)


get_anno_genes = function(res, fwer_threshold=0.05, background=FALSE, go_ids=NULL, genes=NULL){
	
	## Check input	
	# given res as input:
	if(!(missing(res))){
		## checks
		# check that it is really a res-object
		if(!(is.list(res)) || is.null(names(res)) || !(all(names(res) == c("results","genes")))){
			stop("Please use an object returned from go_enrich as input (list with 2 elements).\n Alternatively go_ids need to be defined.")
		}
		if(!is.null(go_ids)){
			warning("Unused argument: structure_ids")
		}		
		if(!is.null(genes)){
			warning("Unused argument: genes")
		}
		## get significant GOs given the fwer_threshold
		message("extract enriched GOs...")
		fwers = res[[1]]
		go_ids = fwers[fwers[,7] < fwer_threshold, "node_id"]
		if(length(go_ids) == 0){
			warning(paste("No significantly enriched GOs at FWER-threshold ",fwer_threshold,sep=""))
			return(NULL)
		}
		## get genes
		res_genes = res[[2]]
	# without res as input	
	} else {
		## checks
		# GO-IDs
		if(is.null(go_ids)){
			stop("Please define go_ids.")
		}
		# background
		if(background){
			warning("Unused argument: background")
		}
	}

	## GO-graph and GO-annotations	
	# find children
	message("find child nodes of GOs...")
	children = list()
	for(go in go_ids){
		# go_id in ontology
		go_id_term = term[term[,4]==go, 1]
		if(length(go_id_term) == 0){ 
			next
		}	
		# go-categories that are children of *go* (cols graph_path: id, parent, child, relation, distance, ...)
		go_children = graph_path[graph_path[,2]==go_id_term,]
		# GO-IDs of child nodes (and the node itself)
		child_gos = term[term[,1] %in% go_children[,3],4]
		children[[go]] = child_gos
	}
	
	# 3) find genes annotated to child nodes
	message("find genes annotated to child nodes...")
	# remove obsolete terms
	term = term[term[,5]==0,]
	# remove empty genes in annotations
	go_anno = go_anno[go_anno[,2]!="",]
	# restrict annotations to GOs present in term
	anno = go_anno[go_anno[,1] %in% term[,4],]
	# restrict to genes
	if (!missing(res)){
		if(!background && all(res_genes %in% c(0,1))){ # reduce to candidate genes (hyper, background=FALSE)
			genes = names(res_genes[res_genes==1])
		} else if(background & all(res_genes==1)){ # dont reduce (hyper, backgr=TRUE, no bg defined)
			genes = NULL
		} else { # reduce to input genes (hyper, background=TRUE and defined, wilcoxon)
			genes = names(res_genes)
		}
	}
	if(!is.null(genes)){
		anno = anno[anno[,2] %in% genes, ]  
	}
	if(nrow(anno)==0){
		warning("No GO-annotations found for the input genes.")
		return(NULL)
	}
	# restrict to child-GOs
	child_anno = lapply(children, function(x) unique(anno[anno[,1] %in% x,2]))
	# rearrange to dataframe
	message("rearrange output...")
	go_id = lapply(seq_along(child_anno), function(i) rep(names(child_anno)[i], length(child_anno[[i]])))
	go_id = unlist(go_id)
	if(length(go_id)==0){
		warning("No GO-annotations found for the input.")
		return(NULL)
	}
	anno_gene = unlist(child_anno)
	out = data.frame(go_id, anno_gene)
	
	if(!(missing(res))){
		fwer = fwers[match(out[,1], fwers[,2]),7]
		score = res_genes[as.character(out$anno_gene)]
		# replace NA with 0 for background genes
		score[is.na(score)] = 0
		out = cbind(out, fwer, score)		
		out = out[order(out$FWER, out$go_id, out$score, out$anno_gene),]
	} else {
		out = out[order(out$go_id, out$anno_gene),]
	}
		
	row.names(out) = 1:nrow(out)
	message("Done.")

	return(out)
}	
