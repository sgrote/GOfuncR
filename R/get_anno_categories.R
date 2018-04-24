

# return a dataframe with genes annotated to GOs 

# input: genes, optional: database, alternative: annotation dataframe;
#     optional term-object or term-table directory

# output: gene, go_id


get_anno_categories = function(genes, database="Homo.sapiens", annotations=NULL, term_df=NULL, godir=NULL, silent=FALSE){
    
    if (!missing(genes)){
        genes = as.character(genes)
    }
    
    ### get annotations
    
    if (is.null(annotations)){
        
        ## A) annotation package
        load_db(database, silent)
        # if genes are not provided use all from database (useful for default background in hypergeometric test)
        if (missing(genes)){
            genes = keys(get(database), keytype="SYMBOL")
        } else if(!any(genes %in% keys(get(database), keytype="SYMBOL"))){
            stop("None of the genes entered are present in the SYMBOL column of '" ,database, "'. Check head(keys(", database, ", keytype='SYMBOL')) to see valid examples.")
        }
        
        ## find annotated GO-categories
        if (!silent){
            message("Find associated categories using database '",database,"'...")
        }
        # load GO-annotation
        go_anno = suppressMessages(select(get(database), keys=genes, columns=c("SYMBOL","GO"), keytype="SYMBOL"))
        go_anno = go_anno[!is.na(go_anno[,2]), 1:2]
    
    } else {
        
        ## B) custom annotation dataframe
        if (!silent){
            message("Find associated categories using custom annotations...")
        }      
        if (missing(genes)){
            go_anno = annotations
        } else {
            go_anno = annotations[annotations[,1] %in% genes, ]
        }
    }

    ### reduce to go-graph
    
    if (!silent){
        message("Remove annotated categories not present in GO-graph...")
    }      
    
    # check if term is user-defined or get the integrated version
    term = eval_term(term_df, godir)
    
    # remove obsolete terms
    term = term[term[,5]==0,]
    # restrict annotations to GOs present in term
    out = go_anno[go_anno[,2] %in% term[,4],]

    if(nrow(out)==0){
        stop("No valid GO-annotations found for the input genes.")
    }
    # create output
    colnames(out) = c("gene", "go_id")
    out = unique(out)

    # sort
    out = out[order(out$gene, out$go_id),]
    rownames(out) = 1:nrow(out)

    return(out)
}   


