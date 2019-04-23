

## for hyper

refine = function(res, pval){

    # TODO: remove hardcoded raw_p_overrep

    # get significant GO-categories
    stats = res[[1]]
    signi_stats = stats[stats$raw_p_overrep < pval, ]
    message("Found ", nrow(signi_stats), " significant categories at p-value threshold ", pval)

    # infer ontology
    # TODO remove redundancy with get_anno_scores
    onto = res[[3]][res[[3]][,1] == "go_graph", ]
    if (onto[1,2] == "custom"){
        godir = onto[1,3]
        term = read.table(paste0(godir,"/term.txt"),sep="\t",quote="",comment.char="",as.is=TRUE)
        graph_path = read.table(paste0(godir,"/graph_path.txt"),sep="\t",quote="",comment.char="",as.is=TRUE)
    }

    # get scores for nodes
    anno = get_anno_scores(res, signi_stats$node_id, annotations=NULL, go_roots_only=TRUE)

    # add term-ID
    anno = cbind(anno, term_id=go_to_term(anno$go_id, term))

    # add GO-ID of children
    graph_path$go_id = term_to_go(graph_path[,3], term)
    # intialize refinement results
    refined = data.frame(node_id=signi_stats$node_id, raw_p=signi_stats$raw_p_overrep, refined_p=NA) 
    # loop over roots
    # TODO: for contingency root nodes are not part of anno, have to get them differently
    root_ids = unique(anno$root_id)
    for (r in root_ids){
        message("\nRefine root node ", r, "...")
        # restrict graph to significant+root categories from this root node
        anno_one_root = anno[anno$root_id == r,]
        term_ids = unique(anno_one_root$term_id)
        sub_graph_path = graph_path[graph_path[,2] %in% term_ids & graph_path[,3] %in% term_ids, ]
        # remove dist=0
        sub_graph_path = sub_graph_path[sub_graph_path[,5] != 0,]
        # divide into root and significant-node annoatios
        anno_signi = anno_one_root[anno_one_root$go_id != r,]
        anno_root = anno_one_root[anno_one_root$go_id == r,]
        # aggregate scores per root (stable across refinement rounds)
        scores_root = tapply(anno_root[,4], anno_root[,1], function(x) c(sum(x[]), length(x)-sum(x)))
        scores_root = do.call(rbind, scores_root)
        # recursively compute refinement (update 'refined')
        refined = refine_algo(anno_signi, scores_root, sub_graph_path, pval, refined)
    }

    # merge with original data, all=T just to be sure, should always have the same
    out = merge(signi_stats, refined, by="node_id", all=T)
    out$signif = out$refined_p < pval

    message("\nDone.")

    return(out)
}


# recursive algorithm for refinement
refine_algo = function(anno_signi, scores_root, sub_graph_path, pval, refined){
    
    # get leaves sub_graph_path (id, parent, child, ...)
    leaves = unique(sub_graph_path[!(sub_graph_path[,3] %in% sub_graph_path[,2]), "go_id"])
    message("Found ", length(leaves), " leaves:")
    print(leaves)
    
    # return updated p-vals table if no leaves are left
    # (also works if dist=0 is not removed from graph_path)
    if (length(leaves)==0){
        return(refined)
    }
    
    # aggregate annotations for leaves
    anno_leaves = anno_signi[anno_signi$go_id %in% leaves, ]
    empty_leaves = leaves[!leaves %in% anno_leaves[,1]]
    names(empty_leaves) = empty_leaves

    # TODO: from here it's hyper-specific

    # counts of 1 and 0 genes in a node
    scores_leaves = tapply(anno_leaves[,4], anno_leaves[,1], function(x) c(sum(x[]), length(x)-sum(x)))
    # add leaves without any genes left - if there's none doesn't matter
    scores_leaves = c(scores_leaves, lapply(empty_leaves, function(x) c(0,0)))
    # collapse
    scores_leaves = data.frame(go_id=names(scores_leaves), do.call(rbind, scores_leaves))
    
    # save time: run test for each combination of scores only once
    scores_leaves$score_id = paste(scores_leaves[,2], scores_leaves[,3], sep="_")
    unique_scores = unique(scores_leaves[,2:4])
    # hyper(candi_node, candi_root, bg_node, bg_root, under=FALSE){
    new_p = hyper(unique_scores[,1], scores_root[,1], unique_scores[,2], scores_root[,2])

    # TODO: from here it's general again

    # add to leaves
    scores_leaves$new_p = new_p[match(scores_leaves$score_id, unique_scores$score_id)]
    
    # add to output
    refined[match(scores_leaves$go_id, refined$node_id), "refined_p"] = scores_leaves$new_p

    # extract leaves that are still significant, if there are none it doesn't matter
    signi_leaves = scores_leaves[scores_leaves$new_p < pval, ]
    #signi_leaves$term_id = anno_signi[match(signi_leaves$go_id, anno_signi$go_id), "term_id"]
    message(nrow(signi_leaves), " of the ", length(leaves), " leaves stay significant after refinement.")
    
    # remove their genes from all their ancestors
    #if (nrow(signi_leaves) > 1) message("node-gene annotations:")
    for (i in seq_along(signi_leaves[,1])){
        sl = signi_leaves[i,]
        genes = anno_signi[anno_signi$go_id == sl$go_id, "gene"]
        # term-IDs of all ancestors
        parents = sub_graph_path[sub_graph_path$go_id == sl$go_id, 2]
        to_remove = anno_signi$term_id %in% parents & anno_signi$gene %in% genes
        anno_signi = anno_signi[!to_remove,]
        #print(nrow(anno_signi))
    }
    
    # remove all leaves from graph_path
    sub_graph_path = sub_graph_path[!(sub_graph_path$go_id %in% leaves), ]
    
    # refine next level with updated annotations, graph_path and final 'refined'-table
    message("Refine next level...")
    refine_algo(anno_signi, scores_root, sub_graph_path, pval, refined)
}



