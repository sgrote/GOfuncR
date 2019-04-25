# take a result from go_enrich() and a vector of GO-IDs to plot annotated scores
# plotting type depends on the performed test in go_enrich which is automatically recognized
# (fwer_threshold is not supported anymore, to not discriminate high_A or high_B in binomial, contingency)


plot_anno_scores = function(res, go_ids, annotations=NULL){
    

    ### check input
    check_res(res)
    # check that all go_ids are in res
    if (!all(go_ids %in% res[[1]][,2])){
        inval = go_ids[!go_ids %in% res[[1]][,2]]
        stop("go_ids not present in res (go_enrich result): ", paste(inval,collapse=", "))
    }
    # infer test
    in_genes = res[[2]]
    test = infer_test(in_genes)
    
    # infer ontology
    onto = load_onto(res[[3]])
    term = onto[[1]]
    graph_path = onto[[2]]

    # infer root_nodes
    root_names = rev(sort(unique(res[[1]][,1])))
    # if some (default) root node is skipped in go_enrich, add it for stable colors
    default_roots = c("molecular_function","cellular_component","biological_process")
    if (all(root_names %in% default_roots)){
        root_names = default_roots
    }
    
    # colors and IDs for root nodes
    root_names_id = term[match(root_names, term[,2]) ,]
    pie_cols = c("#F15A60","#7BC36A","#599BD3","#F9A75B","#9E67AB","#CE7058","#D77FB4")
    root_cols = data.frame(root=root_names_id[,4], col=pie_cols[1:nrow(root_names_id)], stringsAsFactors=FALSE)
    
    # reduce to used domains
    # (to not plot root-node for unused domain, still have before for stable colors)
    root_names_id = root_names_id[root_names_id[,2] %in% res[[1]][,1], ]
    root_ids = root_names_id[,4]
    
    # just in case
    go_ids = as.character(go_ids)
    # keep order of input GO's (which gets messed up in get_anno_genes by *apply)
    ordere = data.frame(go_ids, rank=1:length(go_ids))
    
    # get annotations and scores for GO-IDs
    anno_scores = get_anno_scores(res, go_ids, term, graph_path, annotations)

    # aggregate scores in nodes (wilcox: plot score distribution)
    if (test == "hyper"){
        # counts of 1 and 0 genes in a node
        anno_scores = tapply(anno_scores[,3], anno_scores[,1], function(x) c(sum(x), length(x)-sum(x)))
        anno_scores = data.frame(go_id = names(anno_scores), do.call(rbind, anno_scores))
    } else if (test %in% c("binomial", "contingency")){
        # sums of scores in a node (binom + conti)
        anno_scores = aggregate(anno_scores[,3:ncol(anno_scores)], list(go_id=anno_scores[,1]), sum)
    }

    ### get annotation for root nodes   (conti independent of root nodes)
    if (test != "contingency"){
        
        # get annotation and scores for root nodes
        root_anno_scores = get_anno_scores(res, root_ids, term, graph_path, annotations)
        # order alphabetically (which is rev(order(IDs)) for default GO)
        root_anno_scores = root_anno_scores[rev(order(root_anno_scores[,1])), ]

        # aggregate scores in root nodes
        if (test != "wilcoxon"){
            if (test == "hyper"){ 
                # counts of 1 and 0 genes in a node
                root_anno_scores = tapply(root_anno_scores[,3], root_anno_scores[,1], function(x) c(sum(x[]), length(x)-sum(x)))
                root_anno_scores = data.frame(go_id = names(root_anno_scores), do.call(rbind, root_anno_scores))
            } else { 
                # sums of scores in a node (binom + conti)
                root_anno_scores = aggregate(root_anno_scores[,3:ncol(root_anno_scores)], list(go_id=root_anno_scores[,1]), sum)
            }
            # add colors and root_node_name
            root_anno_scores$root_name = root_names_id[match(root_anno_scores[,1], root_names_id[,4]), 2]
            root_anno_scores$root_col = root_cols[match(root_anno_scores[,1], root_cols[,1]), 2]
            # merge nodes with root node info
            matched_root_name = get_names(anno_scores[,1], term)[,3] # get names also has root name
            anno_scores$root_id = root_names_id[match(matched_root_name, root_names_id[,2]), 4]
            anno_scores = cbind(anno_scores, root_anno_scores[match(anno_scores$root_id, root_anno_scores[,1]), 2:ncol(root_anno_scores)])
        } else { 
            # for wilcox leave unaggregated version but create table with median, name, col
            root_info = aggregate(root_anno_scores[,3], list(go_id=root_anno_scores[,1]), median)
            root_info$root_name = get_names(root_info[,1], term)[,2]
            root_info$root_col = root_cols[match(root_info[,1], root_cols[,1]), 2]
        }
    }
    
    # recover original order (aggregate and get_anno_genes sorts output alphabetically)
    anno_scores = anno_scores[order(ordere[match(anno_scores$go_id, ordere$go_ids), 2]),]
    rownames(anno_scores) = 1:nrow(anno_scores)

    # plot and get stats returned
    if (test == "hyper"){
        stats = plot_hyper(anno_scores, root_anno_scores)
    } else if (test == "binomial"){
        stats = plot_binomial(anno_scores, root_anno_scores)
    } else if (test == "contingency"){
        stats = plot_conti(anno_scores)
    } else if (test == "wilcoxon"){
        stats = plot_wilcox(anno_scores, root_anno_scores, root_info, term)
    }
    
    return(invisible(stats))
}
    
    
    

