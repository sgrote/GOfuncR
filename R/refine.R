

# annotations = custom annotations dataframe if needed
refine = function(res, pval, pcol=5, annotations=NULL){
    
    # check input
    check_res(res)
    onto = load_onto(res[[3]])
    term = onto[[1]]
    graph_path = onto[[2]]
    test = infer_test(res[[2]])

    # check pval-column, allow string or number
    if (! is.numeric(pcol)){
        pcol = which(colnames(res[[1]]) == pcol)
    }
    if (length(pcol) == 0 || ! pcol %in% c(4,5)){
        colnam = paste(colnames(res[[1]])[4:5], collapse = "' or '")
        stop("Invalid 'pcol' argument.\nPlease use 4 or 5 or the corresponding colnames '", colnam, "'.")
    }
    pcol_string = colnames(res[[1]])[pcol]
    message("\nRefinement on '", pcol_string, "' using p-value threshold ", pval, "...\n")

    # define direction for refinement test
    # col 4 is underrep/low-scores and col 5 overrep/high-scores/high-A/B in each test
    if (pcol == 4){
        low = TRUE
    } else if (pcol == 5){
        low = FALSE
    }

    # get significant GO-categories
    stats = res[[1]]
    signi_stats = stats[stats[,pcol] < pval, ]
    if (nrow(signi_stats) == 0){
        stop("There are no significant categories at p-value threshold ", pval)
    }
    message("Found ", nrow(signi_stats), " significant categories at p-value threshold ", pval)
    go_ids = signi_stats$node_id
    
    # add root nodes for significant nodes
    if (test != "contingency"){
        root_names = unique(signi_stats[,1])
        root_names_id = term[match(root_names, term[,2]),]
        root_ids = root_names_id[,4]
        go_ids = c(go_ids, root_ids)
        # add root_ids to GO-ids for combining them with output
        matched_root_name = get_names(go_ids, term) # get names also returns root-node-name
        matched_root_name$root_id = root_names_id[match(matched_root_name[,3], root_names_id[,2]), 4]
    }
    
    # get annotated genes with scores
    message("\nGet annotations for categories:")
    anno = get_anno_scores(res, go_ids, term, graph_path, annotations)

    # add term-ID and root_ids
    anno = cbind(anno, term_id=go_to_term(anno$go_id, term))
    if (test != "contingency"){
        anno = cbind(anno, root_id=matched_root_name[match(anno[,1], matched_root_name[,1]), "root_id"])
        anno$root_id = as.character(anno$root_id)
    }

    # add GO-ID of children
    graph_path$go_id = term_to_go(graph_path[,3], term)
    # intialize refinement results
    refined = data.frame(node_id=signi_stats$node_id, raw_p=signi_stats[,pcol], refined_p=NA) 
    # loop over roots
    # TODO: for contingency root nodes are not part of anno, have to get them differently
    for (r in root_ids){
        message("\nRefine root node ", r, "...")
        # restrict graph to significant+root categories from this root node
        anno_one_root = anno[anno$root_id == r,]
        term_ids = unique(anno_one_root$term_id)
        sub_graph_path = graph_path[graph_path[,2] %in% term_ids & graph_path[,3] %in% term_ids, ]
        # remove dist=0
        sub_graph_path = sub_graph_path[sub_graph_path[,5] != 0,]
        # divide into root and significant-node annotations
        anno_signi = anno_one_root[anno_one_root$go_id != r,]
        anno_root = anno_one_root[anno_one_root$go_id == r,]
        # aggregate scores per root (stable across refinement rounds)
        # (leave individual scores for wilcox)
        if (test == "hyper"){
            # counts of 1 and 0 genes in a node
            scores_root = c(sum(anno_root[,3]), length(anno_root[,3])-sum(anno_root[,3]))
        } else if (test == "wilcoxon"){
            # genes, scores
            scores_root = anno_root[,2:3]
        } else if (test == "binomial"){
            scores_root = colSums(anno_root[,3:4])
        } else if (test == "contingency"){
            scores_root = colSums(anno_root[,3:6])
        }
        # recursively compute refinement (update 'refined')
        refined = refine_algo(anno_signi, scores_root, sub_graph_path, pval, refined, test, low, first=TRUE)
    }

    # merge with original data, all=T just to be sure, should always have the same
    out = merge(signi_stats, refined, by="node_id", all=TRUE, sort=FALSE)
    out$signif = out$refined_p < pval
    # TODO remove raw_p column after testing
    colnames(out)[colnames(out) == "refined_p"] = paste0("refined_", substring(pcol_string, 5))
    
    message("\nDone.")

    return(out)
}


# recursive algorithm for refinement, bottom-up
# anno_signi: data.frame(go_id, gene, scores, ..., term_id, ...)
# scores_root: aggregated scores per root, e.g. c(candi, bg) for hyper
# sub_graph_path: data.frame(id, term_id_ancestor, term_id_child, ..., go_id=go_id of child)
#   for all significant nodes, all distances>0
# pval: p-value threshold for significance
# refined: data.frame(node_id, raw_p, refined_p) - output, recursively updated
# test: e.g. "wilcoxon"
# low:  T/F direction of test, e.g. under-/overrep
# first: T/F first round of refinement, i.e. the very leaves of the significant sub-graph
refine_algo = function(anno_signi, scores_root, sub_graph_path, pval, refined, test, low, first){
    
    # get go_id of leaves from sub_graph_path (id, parent, child, ...)
    leaves = unique(sub_graph_path[!(sub_graph_path[,3] %in% sub_graph_path[,2]), "go_id"])
    message("Found ", length(leaves), " leaves:")
    print(leaves)
    
    # return updated p-vals table if no leaves are left
    # (also works if dist=0 is not removed from graph_path)
    if (length(leaves)==0){
        return(refined)
    }
    
    # annotations for leaves
    anno_leaves = anno_signi[anno_signi$go_id %in% leaves, ]
    empty_leaves = leaves[!leaves %in% anno_leaves[,1]]
    
    # run category test for leaves; data.frame(go_id, new_p)
    # TODO: maybe skip in first round after testing that raw_p == leaf_p for all absolute leaves
    if (test == "hyper"){
        new_p_leaves = hyper_nodes(anno_leaves, empty_leaves, scores_root, low)
    } else if (test == "wilcoxon"){
        new_p_leaves = wilcox_nodes(anno_leaves, empty_leaves, scores_root, low)
    } else if (test == "binomial"){
        new_p_leaves = binom_nodes(anno_leaves, empty_leaves, scores_root, low)
    }
    
    # add to output
    refined[match(new_p_leaves$go_id, refined$node_id), "refined_p"] = new_p_leaves$new_p

    # check that p-value for leaves is the same before and after refinement
    if (first){
        abs_leaves = refined[!is.na(refined$new_p),]
        if (! all.equal(abs_leaves$raw_p, abs_leaves$refined_p)){
            print(abs_leaves)
            stop("raw_p != refined_p for lowest level of significant nodes.")
        }
    }

    # extract leaves that are still significant, if there are none it doesn't matter
    signi_leaves = new_p_leaves[new_p_leaves$new_p < pval, ]
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
    refine_algo(anno_signi, scores_root, sub_graph_path, pval, refined, test, low, first=FALSE)
}



