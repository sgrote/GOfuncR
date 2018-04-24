


# return a dataframe with child and parent nodes and distances 
get_parent_nodes = function(go_ids, term_df=NULL, graph_path_df=NULL, godir=NULL){

    # get term and grpah_path depending on input
    onto = eval_onto_input(term_df, graph_path_df, godir)
    term = onto[[1]]
    graph_path = onto[[2]]
    # get term.txt IDs for GO-IDs
    go_ids_term = go_to_term(go_ids, term)
    # go-categories that are parents of *go* (cols graph_path: id, parent, child, relation, distance, ...)
    # term-IDs and distance
    go_parent_id = graph_path[graph_path[,3] %in% go_ids_term,]
    # GO-IDS of graph_path term-IDs
    parent = term_to_go(go_parent_id[,2], term)
    children = term_to_go(go_parent_id[,3], term)
    # add name to parent
    names = get_names(parent, term)[,"go_name"]
    out = data.frame(children, parent, names, go_parent_id[,5], stringsAsFactors=FALSE)
    # remove 'all'-root node
    out = out[out[,2] != "all",]
    out = sort_out(out)
    colnames(out) = c("child_go_id","parent_go_id","parent_name","distance")
    
    return(out)
}



# return a dataframe with parent and child nodes and distances 
get_child_nodes = function(go_ids, term_df=NULL, graph_path_df=NULL, godir=NULL){

    # get term and grpah_path depending on input
    onto = eval_onto_input(term_df, graph_path_df, godir)
    term = onto[[1]]
    graph_path = onto[[2]]
    # get term.txt IDs for GO-IDs
    go_ids_term = go_to_term(go_ids, term)
    # go-categories that are children of *go* (cols graph_path: id, parent, child, relation, distance, ...)
    # term-IDs and distance
    go_children_id = graph_path[graph_path[,2] %in% go_ids_term,]
    # GO-IDS of graph_path term-IDs
    parent = term_to_go(go_children_id[,2], term)
    children = term_to_go(go_children_id[,3], term)
    # add name to child
    names = get_names(children, term)[,"go_name"]
    out = data.frame(parent, children, names, go_children_id[,5], stringsAsFactors=FALSE)
    out = sort_out(out)
    colnames(out) = c("parent_go_id","child_go_id","child_name","distance")

    return(out)
}

# get term.txt IDs for GO-IDs
go_to_term = function(go_ids, term){
    go_ids = as.character(go_ids)
    # remove obsolete terms
    term = term[term[,5]==0,]
    # get term-IDs of GOs
    go_ids_term = term[term[,4] %in% go_ids, 1]
    if(length(go_ids_term) == 0){
        stop("GO-IDs not found in ontology or obsolete.")
    }
    return(go_ids_term)
}

# get GO-IDs of term.txt IDs
term_to_go = function(term_ids, term){
    go_ids = term[match(term_ids, term[,1]) ,4]
    return(go_ids)
}

# rearrange output
sort_out = function(out){
    # sort
    out = out[order(out[,1],out[,4],out[,2],out[,3]),]
    # remove duplicates (ordered by distance, so shortest dist is kept)
    out = out[!(duplicated(paste(out[,1],out[,2]))),]
    rownames(out) = 1:nrow(out)
    return(out)
}
