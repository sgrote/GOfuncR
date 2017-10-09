

# return a dataframe with genes annotated to GOs 

# input: genes, optional: database

# output: gene, go_id


get_anno_categories = function(genes, database="Homo.sapiens"){
    
    ## Check input  
    if(is.null(genes)){
        stop("Please provide gene-symbols as input.")
    }
    if (!suppressPackageStartupMessages(suppressWarnings(require(database, character.only=TRUE)))){
		stop(paste0("database '" ,database, "' is not installed. Please install it from bioconductor."))
	}
    if(!any(genes %in% keys(get(database), keytype="SYMBOL"))){
		stop(paste0("None of the genes entered are present in the SYMBOL column of '" ,database, "'. Check head(keys(", database, ", keytype='SYMBOL')) to see valid examples."))
	}
	
    ## find annotated GO-categories
    message(paste("find associated categories using database '",database,"'...",sep=""))
    # load GO-annotation
    go_anno = select(get(database), keys=genes, columns=c("SYMBOL","GO"), keytype="SYMBOL")
	go_anno = go_anno[!is.na(go_anno[,2]), 1:2]

    ## remove obsolete terms
    term = term[term[,5]==0,]
    # restrict annotations to GOs present in term
    out = go_anno[go_anno[,2] %in% term[,4],]

    if(nrow(out)==0){
        stop("No valid GO-annotations found for the input genes.")
    }
    # create output
    colnames(out) = c("gene", "go_id")

    # add category name 
    # TODO: don't add name here but show get_names example in vignette and man
#    out = cbind(out, get_names(out$go_id)[,c("go_name","root_node")])
#    # sort
    out = out[order(out$go_id, out$gene),]
    rownames(out) = 1:nrow(out)

    return(out)
}   


