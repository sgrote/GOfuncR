# return a dataframe with parent and child nodes and distances 

get_child_nodes = function(go_ids){

    go_ids = as.character(go_ids)
    # remove obsolete terms
    term = term[term[,5]==0,]
    # get term-IDs of GOs
    go_ids_term = term[term[,4] %in% go_ids, 1]
    if(length(go_ids_term) == 0){
        stop("GO-IDs not found in ontology or obsolete.")
    }
    # go-categories that are children of *go* (cols graph_path: id, parent, child, relation, distance, ...)
    # term-IDs and distance
    go_children_id = graph_path[graph_path[,2] %in% go_ids_term,]
    # names and GO-IDS
    parent = term[match(go_children_id[,2],term[,1]) ,4]
    # GO-IDs of child nodes (and the node itself)
    children = term[match(go_children_id[,3],term[,1]) ,4]
    # NEW: add name to child
    names = get_names(children)[,"go_name"]
    out = data.frame(parent, children, names, go_children_id[,5])
    colnames(out) = c("parent_go_id","child_go_id","child_name","distance")
    out[,1:3] = apply(out[,1:3], 2, as.character)
    # sort
    out = out[order(out[,1],out[,4],out[,2],out[,3]),]
    # NEW: remove duplicates (ordered by distance, so shortest dist is kept)
    out = out[!(duplicated(paste(out[,1],out[,2]))),]

    rownames(out) = 1:nrow(out)

    return(out)
}
