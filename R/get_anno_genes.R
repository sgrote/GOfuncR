
# return a dataframe with genes annotated to GOs 

# input: go_ids, optional: database, genes, custom ontology

# output: go_id, gene


get_anno_genes = function(go_ids, database="Homo.sapiens", genes=NULL, annotations=NULL, term_df=NULL, graph_path_df=NULL, godir=NULL){
    
    ## Check input  
    go_ids = as.character(go_ids)
    if (!is.vector(go_ids) || !all(substr(go_ids,1,3) == "GO:")){
        stop("Please provide GO-IDs as input, e.g. go_ids=c('GO:0072221','GO:0004945')")
    }
    
    # get term and grpah_path depending on input
    onto = eval_onto_input(term_df, graph_path_df, godir)
    term = onto[[1]]
    graph_path = onto[[2]]
    
    ### get annotations
    
    # child nodes    
    message("find child nodes of GOs...")
    children = get_child_nodes(go_ids, term, graph_path)[,1:2]
    uni_child = unique(children[,2])
    
    if (is.null(annotations)){
        
        ## A) annotation package
        load_db(database)
        # find genes annotated to child nodes
        message("find genes annotated to child nodes using database '", database,"'...")
        child_anno = suppressMessages(select(get(database), keys=uni_child, columns=c("GO","SYMBOL"), keytype="GO"))
        child_anno = child_anno[!is.na(child_anno[,4]), c(1,4)]
    } else {
        ## B) custom annotation dataframe
        message("Find associated categories using custom annotations...")
        child_anno = annotations[annotations[,2] %in% uni_child, 2:1]
        colnames(child_anno) = c("GO","SYMBOL")
    }
    
    ### reduce to go-graph
    
    # remove obsolete terms
    message("Remove annotated categories not present in GO-graph...")
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


