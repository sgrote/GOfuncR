

#### hypergeometric test

hyper = function(candi_node, candi_root, bg_node, bg_root, under=FALSE){
    total_node = candi_node + bg_node
    if (under){
        p = phyper(candi_node, candi_root, bg_root, total_node)
    } else {
        p = phyper(candi_node-1, candi_root, bg_root, total_node, lower.tail=FALSE)
    }
    return(p)
}

# to test all leaves at once in refinement
# also hande empty_nodes = nodes with no annotations left during refinement steps 
# anno_nodes: data.frame(go_id, gene, score 0/1)
# empty_nodes: vector of go-ids of empty nodes
# scores_root: data.frame(go_id, scores)
# low: T/F under-/overrep
# out: go_id, new_p
hyper_nodes = function(anno_nodes, empty_nodes, scores_root, low=FALSE){
    
    # counts of 1 and 0 genes in a node
    scores_nodes = tapply(anno_nodes[,3], anno_nodes[,1], function(x) c(sum(x[]), length(x)-sum(x)))
    # add nodes without any genes left - if there's none doesn't matter
    names(empty_nodes) = empty_nodes
    scores_nodes = c(scores_nodes, lapply(empty_nodes, function(x) c(0,0)))
    # collapse
    scores_nodes = data.frame(go_id=names(scores_nodes), do.call(rbind, scores_nodes), row.names=NULL)
    
    # save time: run test for each combination of scores only once
    scores_nodes$score_id = paste(scores_nodes[,2], scores_nodes[,3], sep="_")
    unique_scores = unique(scores_nodes[,2:4])
    # hyper(candi_node, candi_root, bg_node, bg_root, under=FALSE){
    new_p = hyper(unique_scores[,1], scores_root[,1], unique_scores[,2], scores_root[,2], low)
    
    # add to leaves
    scores_nodes$new_p = new_p[match(scores_nodes$score_id, unique_scores$score_id)]
    out = scores_nodes[,c("go_id", "new_p")]
    return(out)
}


#### wilcoxon rank sum test

# input: data.frames [gene, score]
wilcox = function(anno_genes_node, anno_genes_root, low=FALSE){
    scores_outside_node = anno_genes_root[!(anno_genes_root[,1] %in% anno_genes_node[,1]), 2]
    scores_node = anno_genes_node[,2]
    if (low){
        willi = wilcox.test(scores_node, scores_outside_node, alternative="less", exact=FALSE)
    } else {
        willi = wilcox.test(scores_node, scores_outside_node, alternative="greater", exact=FALSE)
    }
    return(willi$p.value)
}



