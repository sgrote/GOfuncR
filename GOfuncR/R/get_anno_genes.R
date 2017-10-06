
# return a dataframe with genes annotated to GOs 

# input: go_ids, optional: genes, ref_genome

# output: go_id, gene


get_anno_genes = function(go_ids, ref_genome="grch37", genes=NULL){
    
    ## Check input  
    # GO-IDs
    if (!is.vector(go_ids) || !all(substr(go_ids,1,3) == "GO:")){
        stop("Please provide GO-IDs as input, e.g. go_ids=c('GO:0072221','GO:0004945')")
    }
    # reference genome
    if (!(ref_genome %in% c("grch37","grch38","grcm38"))){
        stop("Please use ref_genome='grch37', ref_genome='grch38' or ref_genome='grcm38'")
    }

    ## GO-graph and GO-annotations  
    # find children
    message("find child nodes of GOs...")
    children = get_child_nodes(go_ids)
    children = tapply(children[,2], children[,1], as.character, simplify=FALSE)

    # 3) find genes annotated to child nodes
    message(paste("find genes annotated to child nodes using ref_genome ",ref_genome,"...",sep=""))
    # remove obsolete terms
    term = term[term[,5]==0,]
    # load GO-annotation
    go_anno = get(paste("go_anno_", ref_genome, sep=""))
    # remove empty genes in annotations
    go_anno = go_anno[go_anno[,2]!="",]
    # restrict annotations to GOs present in term
    anno = go_anno[go_anno[,1] %in% term[,4],]
    # restrict to genes
    if(!is.null(genes)){
        anno = anno[anno[,2] %in% genes, ]  
    }
    if(nrow(anno)==0){
        warning("No GO-annotations found for the input genes.")
        return(NULL)
    }
    # restrict to child-GOs
    child_anno = lapply(children, function(x) unique(anno[anno[,1] %in% x,2]))
    # rearrange to dataframe
    message("rearrange output...")
    go_id = lapply(seq_along(child_anno), function(i) rep(names(child_anno)[i], length(child_anno[[i]])))
    go_id = unlist(go_id)
    if(length(go_id)==0){
        warning("No GO-annotations found for the input.")
        return(NULL)
    }
    anno_gene = unlist(child_anno)
    out = data.frame(go_id, anno_gene)

    # sort  
    out = out[mixedorder(out$anno_gene),]
    out = out[order(out$go_id),]

    out[,1:2] = apply(out[,1:2], 2, as.character)
    row.names(out) = 1:nrow(out)
    message("Done.")

    return(out)
}   
