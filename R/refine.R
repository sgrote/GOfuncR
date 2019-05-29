

# annotations = custom annotations dataframe if needed
refine = function(res, fwer=0.05, fwer_col=7, annotations=NULL){
    
    # check input
    check_res(res)
    onto = load_onto(res[[3]])
    term = onto[[1]]
    graph_path = onto[[2]]
    test = infer_test(res[[2]])
    if (fwer <= 0 || fwer >= 1){
        stop("Please define a 'fwer' between 0 and 1 (exclusive)")
    }

    # check FWER-column
    if (! is.numeric(fwer_col)){
        fwer_col = which(colnames(res[[1]]) == fwer_col)
    }
    if (length(fwer_col) == 0 || ! fwer_col %in% c(6,7)){
        colnam = paste(colnames(res[[1]])[6:7], collapse = "' or '")
        stop("Invalid 'fwer_col' argument.\nPlease use 6 or 7 or the corresponding colnames '", colnam, "'.")
    }

    # define direction for refinement test
    # col 4 is underrep/low-scores and col 5 overrep/high-scores/high-A/B in each test
    if (fwer_col == 6){
        pcol = 4
        min_pcol = 2
        low = TRUE
    } else if (fwer_col == 7){
        pcol = 5
        min_pcol = 3
        low = FALSE
    }
    pcol_string = colnames(res[[1]])[pcol]

    # get significant GO-categories
    stats = res[[1]]
    signi_stats = stats[stats[,fwer_col] < fwer, ]
    if (nrow(signi_stats) == 0){
        stop("There are no significant categories at FWER threshold ", fwer)
    }
    message("Found ", nrow(signi_stats), " significant categories at FWER threshold ", fwer)
    go_ids = signi_stats$node_id
    
    # get annotated genes with scores for significant GO-categories
    message("\nGet annotations for categories:")
    anno = get_anno_scores(res, go_ids, term, graph_path, annotations)

    # add term-ID
    anno = cbind(anno, term_id=go_to_term(anno$go_id, term))
    
    # add root-IDs
    matched_root_name = get_names(go_ids, term) # get_names also returns root-node-name
    matched_root_name$root_id = term[match(matched_root_name$root_node, term[,2]), 4]
    anno = cbind(anno, root_id=matched_root_name[match(anno[,1], matched_root_name[,1]), "root_id"])
    anno$root_id = as.character(anno$root_id)
    root_ids = unique(matched_root_name$root_id)
    roots = get_names(root_ids, term)
    
    # get annotated genes with scores for root nodes
    if (test != "contingency"){
        message("\nGet annotations for root nodes:")
        anno_root = get_anno_scores(res, root_ids, term, graph_path, annotations)
    }
    # aggregate scores for root nodes
    if (test == "hyper"){
        # counts of 1 and 0 genes in a node
        scores_root_nodes = tapply(anno_root[,3], anno_root[,1],
            function(x){c(sum(x), length(x)-sum(x))}, simplify=FALSE)
    } else if (test == "wilcoxon"){
        # genes, scores
        scores_root_nodes = by(anno_root, anno_root[,1], function(x){x[,c(2,3)]})
    } else if (test == "binomial"){
        scores_root_nodes = by(anno_root, anno_root[,1], function(x){colSums(x[,3:4])})
    } else if (test == "contingency"){
        scores_root_nodes = sapply(root_ids, function(x){NA}, simplify=FALSE)
    }
    
    # add GO-ID of children
    graph_path$go_id = term_to_go(graph_path[,3], term)
    
    # intialize refinement results
    refined = data.frame(node_id=signi_stats$node_id, raw_p=signi_stats[,pcol], refined_p=NA) 
    # loop over roots
    for (i in seq_len(nrow(roots))){
        r = roots[i, 1]
        message("\nRefine root node ", roots[i, 2], "")
        # restrict graph to significant+root categories from this root node
        anno_one_root = anno[anno$root_id == r,]
        # leave root node in sub_graph_path
        term_ids = unique(c(anno_one_root$term_id, go_to_term(r, term)))
        sub_graph_path = graph_path[graph_path[,2] %in% term_ids & graph_path[,3] %in% term_ids, ]
        # remove dist=0
        sub_graph_path = sub_graph_path[sub_graph_path[,5] != 0,]
        # scores per root (stable across refinement rounds)
        scores_root = scores_root_nodes[[r]]
        # convert FWER to p
        min_p = res[[4]][res[[4]][,1] == roots[i, 2], min_pcol]
        pval = fwer_to_p(fwer, min_p)
        message("using p-value threshold ", pval, "...\n")
        # recursively compute refinement (update 'refined')
        refined = refine_algo(anno_one_root, scores_root, sub_graph_path, pval, refined, test, low, first=TRUE)
    }

    # merge with original data, all=T just to be sure, should always have the same
    out = merge(signi_stats, refined, by="node_id", all=TRUE, sort=FALSE)
    out$signif = out$refined_p < pval
    # remove raw_p column after testing
    out = out[, colnames(out) != "raw_p"] 
    colnames(out)[colnames(out) == "refined_p"] = paste0("refined_", substring(pcol_string, 5))
    
    message("\nDone.")

    return(out)
}

fwer_to_p = function(fwer, min_p){
    
    # interpolate with approx
    min_p = sort(min_p)
    unique_p = unique(min_p)
    fwers = ecdf(min_p)(unique_p)
    # avoid NA if fwer too low (just min-p since p_category < p is tested (and not <=))
    if (fwer < min(fwers)){
        p = min(min_p)
    }
    else {
        p = approx(fwers, unique_p, fwer)$y
    }
    return(p)
}


# recursive algorithm for refinement, bottom-up
# anno_signi: data.frame(go_id, gene, scores, ..., term_id, ...)
# scores_root: aggregated scores per root, e.g. c(candi, bg) for hyper,  NA for conti
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
    message("Found ", length(leaves), " leaves...")
    #print(leaves)
    
    # return updated p-vals table if no leaves are left
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
    } else if (test == "contingency"){
        new_p_leaves = conti_nodes(anno_leaves, empty_leaves, low)
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



