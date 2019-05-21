

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
# also handle empty_nodes = nodes with no annotations left during refinement steps 
# anno_nodes: data.frame(go_id, gene, score 0/1, ...)
# empty_nodes: vector of go-ids of empty nodes
# scores_root: c(sum_candi, sum_bg)
# low: T/F under-/overrep
# out: go_id, new_p
hyper_nodes = function(anno_nodes, empty_nodes, scores_root, low=FALSE){
    
    # no annotations at all
    if (nrow(anno_nodes) == 0){
        out = data.frame(go_id=empty_nodes, new_p=1, stringsAsFactors=FALSE)
        return(out)
    }
    # counts of 1 and 0 genes in a node
    scores_nodes = tapply(anno_nodes[,3], anno_nodes[,1], function(x) c(sum(x[]), length(x)-sum(x)))
    # collapse
    scores_nodes = data.frame(go_id=names(scores_nodes), do.call(rbind, scores_nodes),
        stringsAsFactors=FALSE, row.names=NULL)
    
    # save time: run test for each combination of scores only once
    scores_nodes$score_id = paste(scores_nodes[,2], scores_nodes[,3], sep="_")
    unique_scores = unique(scores_nodes[,2:4])
    # hyper(candi_node, candi_root, bg_node, bg_root, under=FALSE){
    new_p = hyper(unique_scores[,1], scores_root[1], unique_scores[,2], scores_root[2], low)
    
    # add p-val to leaves
    scores_nodes$new_p = new_p[match(scores_nodes$score_id, unique_scores$score_id)]
    out = scores_nodes[,c("go_id", "new_p")]
    # add empty leaves
    if (length(empty_nodes) > 0){
        empty_out = data.frame(go_id=empty_nodes, new_p=1)
        out = rbind(out, empty_out)
    }
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


# to test all leaves at once in refinement
# also handle empty_nodes = nodes with no annotations left during refinement steps 
# anno_nodes: data.frame(go_id, gene, score, ...)
# empty_nodes: vector of go-ids of empty nodes
# scores_root: data.frame(gene, score)
# low: T/F lowrank/highrank
# out: go_id, new_p
wilcox_nodes = function(anno_nodes, empty_nodes, scores_root, low=FALSE){

    # no annotations at all
    if (nrow(anno_nodes) == 0){
        out = data.frame(go_id=empty_nodes, new_p=1, stringsAsFactors=FALSE)
        return(out)
    }
    # wilcox test per node
    new_p = by(anno_nodes[,2:3], anno_nodes[,1], wilcox, scores_root, low)
    new_p = do.call(rbind, as.list(new_p))
    # generate output, add empty nodes p=1
    out = data.frame(go_id=rownames(new_p), new_p=new_p[,1], stringsAsFactors=FALSE, row.names=NULL)
    if (length(empty_nodes) > 0){
        empty_out = data.frame(go_id=empty_nodes, new_p=1)
        out = rbind(out, empty_out)
    }
    return(out)
}



#### binomial test

# binom.test(x, n, p)
# x=successes, n=trials, p=prob. of success

binom = function(a_node, b_node, a_root, b_root, low=FALSE){
    p_a = a_root / (a_root + b_root)
    n = a_node + b_node
    if (low){
        # high_B
        bino = binom.test(a_node, n, p=p_a, alternative="less")
    } else {
        # high_A
        bino = binom.test(a_node, n, p=p_a, alternative="greater")
    }
    return(bino$p.value)
}


# to test all leaves at once in refinement
# also handle empty_nodes = nodes with no annotations left during refinement steps
# anno_nodes: data.frame(go_id, gene, counts_A, counts_B, ...)
# empty_nodes: vector of go-ids of empty nodes
# scores_root: c(counts_A, counts_B) in root
# low: T/F high-B/high-A
# out: go_id, new_p
binom_nodes = function(anno_nodes, empty_nodes, scores_root, low=FALSE){

    # no annotations at all
    if (nrow(anno_nodes) == 0){
        out = data.frame(go_id=empty_nodes, new_p=1, stringsAsFactors=FALSE)
        return(out)
    }
    # binomial test per node
    scores_nodes = aggregate(anno_nodes[,3:4], list(go_id=anno_nodes[,1]), sum)
    # save time: run test for each combination of scores only once
    scores_nodes$score_id = paste(scores_nodes[,2], scores_nodes[,3], sep="_")
    unique_scores = unique(scores_nodes[,2:4])
    # binom(a_node, b_node, a_root, b_root, low=FALSE)
    # binom.test doesn't take vectors as input
    new_p = mapply(binom, unique_scores[,1], unique_scores[,2],
        MoreArgs=list(a_root=scores_root[1], b_root=scores_root[2], low=low))
    # add p-val to leaves
    scores_nodes$new_p = new_p[match(scores_nodes$score_id, unique_scores$score_id)]
    out = scores_nodes[,c("go_id", "new_p")]
    # add empty leaves
    if (length(empty_nodes) > 0){
        empty_out = data.frame(go_id=empty_nodes, new_p=1)
        out = rbind(out, empty_out)
    }
    return(out)
}




#### contingency table

# low: T/F p for high_CD/high_AB
conti = function(a_node, b_node, c_node, d_node, low=FALSE){
    # matrix is filled col-wise
    conti_tab = matrix(c(a_node, b_node, c_node, d_node), ncol=2)
    # FUNC original: if one number < 10 --> Fisher's exact test
    if (any(c(a_node, b_node, c_node, d_node) < 10)){
        p = fisher.test(conti_tab)$p.value
    } else {
        p = chisq.test(conti_tab, correct=FALSE)$p.value
    }
    # FUNC original: if A/B > C/D p-val is assigned high_A/B, else high_C/D; other one is set to 1
    if ((a_node / b_node) > (c_node / d_node)){
        p_ab = p
        p_cd = 1
    } else {
        p_ab = 1
        p_cd = p
    }
    if (low){
        return(p_cd)
    } else {
        return(p_ab)
    }
}


# to test all leaves at once in refinement
# also handle empty_nodes = nodes with no annotations left during refinement steps
# anno_nodes: data.frame(go_id, gene, counts_A, counts_B, counts_C, counts_D)
# empty_nodes: vector of go-ids of empty nodes
# low: T/F high C/D  /  high-A/B
# out: go_id, new_p
conti_nodes = function(anno_nodes, empty_nodes, low=FALSE){

    # no annotations at all
    if (nrow(anno_nodes) == 0){
        out = data.frame(go_id=empty_nodes, new_p=1, stringsAsFactors=FALSE)
        return(out)
    }
    # contingency test per node
    scores_nodes = aggregate(anno_nodes[,3:6], list(go_id=anno_nodes[,1]), sum)
    # save time: run test for each combination of scores only once
    scores_nodes$score_id = paste(scores_nodes[,2], scores_nodes[,3], scores_nodes[,4], scores_nodes[,5], sep="_")
    unique_scores = unique(scores_nodes[,2:6])
    # conti(a_node, b_node, c_node, d_node, low=FALSE)
    new_p = mapply(conti, unique_scores[,1], unique_scores[,2], unique_scores[,3], unique_scores[,4],
        MoreArgs=list(low=low))
    # add p-val to leaves
    scores_nodes$new_p = new_p[match(scores_nodes$score_id, unique_scores$score_id)]
    out = scores_nodes[,c("go_id", "new_p")]
    # add empty leaves
    if (length(empty_nodes) > 0){
        empty_out = data.frame(go_id=empty_nodes, new_p=1)
        out = rbind(out, empty_out)
    }
    return(out)
}


