
# return a dataframe with genes annotated to GOs 

# input: go_ids, optional: database, genes

# output: go_id, gene


get_anno_genes = function(go_ids, database="Homo.sapiens", genes=NULL){
    
    ## Check input  
    # GO-IDs
    if (!is.vector(go_ids) || !all(substr(go_ids,1,3) == "GO:")){
        stop("Please provide GO-IDs as input, e.g. go_ids=c('GO:0072221','GO:0004945')")
    }
	# database
	message(paste0("load database '", database, "'..."))
    if (!suppressPackageStartupMessages(suppressMessages(require(database, character.only=TRUE)))){
		stop(paste0("database '" ,database, "' is not installed. Please install it from bioconductor."))
	}

    # child nodes 
    message("find child nodes of GOs...")
    children = get_child_nodes(go_ids)[,1:2]

    # find genes annotated to child nodes
    message(paste("find genes annotated to child nodes using database '", database,"'...",sep=""))
    child_anno = suppressMessages(select(get(database), keys=unique(children[,2]), columns=c("GO","SYMBOL"), keytype="GO"))
    child_anno = child_anno[!is.na(child_anno[,4]), c(1,4)]
    
    # remove obsolete terms
    term = term[term[,5]==0,]
    # restrict annotations to GOs present in term
    child_anno = child_anno[child_anno[,1] %in% term[,4],]
    # restrict to genes
    if(!is.null(genes)){
        child_anno = child_anno[child_anno[,2] %in% genes, ]  
    }
    if(nrow(child_anno) == 0){
        stop("No GO-annotations found for the input genes.")
    }
    # remap to parents (accounts also for multiple parents of same child)
    # (another way would be to get annotations separately for every node, but that might involve recurrent db-queries of shared child-nodes)
	out = merge(children,child_anno, by.x="child_go_id", by.y="GO")
	out = out[,2:3]
	out = unique(out) # remove dupicates (e.g. different children with same annotations)
	colnames(out) = c("go_id", "gene")
	
	
    # sort  
    out = out[order(out$go_id, out$gene),]
    row.names(out) = 1:nrow(out)
    message("Done.")

    return(out)
}   


