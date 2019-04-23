

# hypergeometric test

hyper = function(candi_node, candi_root, bg_node, bg_root, under=FALSE){
    total_node = candi_node + bg_node
    if (under){
        p = phyper(candi_node, candi_root, bg_root, total_node)
    } else {
        p = phyper(candi_node-1, candi_root, bg_root, total_node, lower.tail=FALSE)
    }
    return(p)
}



