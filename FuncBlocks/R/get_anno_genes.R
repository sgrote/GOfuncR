

# return a dataframe with genes annotated to (enriched) GOs 

# input: result list produced by go_enrich, optional: FWER, background
# alternative input: go_ids, optional: genes

# output: go_id, gene, (FWER, candiate/score)


get_anno_genes = function(res, fwer_threshold=0.05, background=FALSE, go_ids=NULL, genes=NULL, ref_genome=NULL){
	
	## Check input	
	# given res as input:
	if(!(missing(res))){
		## checks
		# check that it is really a res-object
		if(!(is.list(res)) || is.null(names(res)) || !(all(names(res) == c("results","genes")))){
			stop("Please use an object returned from go_enrich as input (list with 2 elements).\n Alternatively go_ids need to be defined.")
		}
		if(!is.null(go_ids)){
			warning("Unused argument: go_ids")
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
		## get reference genome
		ref_genome = res[[3]]
	# without res as input	
	} else {
		## checks
		# GO-IDs
		if(is.null(go_ids)){
			stop("Please define go_ids.")
		}
		# reference genome
		if(is.null(ref_genome)){
			stop("Please define a ref_genome ('grch37','grch38','grcm38').")
		}
		# background
		if(background){
			warning("Unused argument: background")
		}
	}
	if (!(ref_genome %in% c("grch37","grch38","grcm38"))){
		stop("Please use 'ref_genome=grch37', 'ref_genome=grch38' or 'ref_genome=grcm38'")
	}

	## GO-graph and GO-annotations	
	# find children
	message("find child nodes of GOs...")
	children = get_child_nodes(go_ids)
	children = tapply(children[,2], children[,1], as.character, simplify=FALSE)
	# uniqueness is not required, takes long. duplicates due to different distances to same child node
#	children = tapply(children[,2], children[,1], function(x){as.character(unique(x))}, simplify=FALSE)


	# 3) find genes annotated to child nodes
	message("find genes annotated to child nodes...")
	# remove obsolete terms
	term = term[term[,5]==0,]
	# load GO-annotation
	go_anno = get(paste("go_anno_", ref_genome, sep=""))
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
		out = out[order(out$fwer, out$go_id, out$score, out$anno_gene),]
	} else {
		out = out[order(out$go_id, out$anno_gene),]
	}
		
	row.names(out) = 1:nrow(out)
	message("Done.")

	return(out)
}	
