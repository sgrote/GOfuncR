
# slimmed clone of plot_anno_scores to get unaggregated annotated genes with scores

# TODO: maybe remove some of the redundancy, use output of this in plot_anno_scores
#       maybe just make helper-functions from bits of it and use them in both plot_anno_scores
#       and refinement, remove this whole file here

get_anno_scores = function(res, go_ids, annotations=NULL, go_roots_only=TRUE){
    
    ### check input
    # check that res could be go_enrich-output
    # TODO: make this a helper-function
    if (!(is.list(res) && all(names(res) == c("results","genes","databases")))){
        stop("Please use an object returned from go_enrich as input (list with 3 elements).")
    }
    # check that all go_ids are in res
    if (!all(go_ids %in% res[[1]][,2])){
        inval = go_ids[!go_ids %in% res[[1]][,2]]
        stop("go_ids not present in res (go_enrich result): ", paste(inval,collapse=", "))
    }
    # infer test
    # TODO: make this a helper-function
    in_genes = res[[2]]
    if (ncol(in_genes) == 2){
        if(all(in_genes[,2] %in% c(1,0))){
            test = "hyper"
        } else {
            test = "wilcoxon"
        }
    } else if(ncol(in_genes) == 3){
        test = "binomial"
    } else if(ncol(in_genes) == 5){
        test = "contingency"
    } else {
        stop("Identification of test failed.")
    }
    
    # infer ontology
    onto = res[[3]][res[[3]][,1] == "go_graph", ]
    if (onto[1,2] == "custom"){
        godir = onto[1,3]
        term = read.table(paste0(godir,"/term.txt"),sep="\t",quote="",comment.char="",as.is=TRUE)
        graph_path = read.table(paste0(godir,"/graph_path.txt"),sep="\t",quote="",comment.char="",as.is=TRUE)
        root_names = rev(sort(unique(res[[1]][,1])))
    } else {
        # use fixed root names for stable colors (if some root node is skipped in go_enrich)
        root_names = c("molecular_function","cellular_component","biological_process")
    }

    # reduce to used domains (either of go_ids or res in general)
    # (to not plot root-node for unused domain, still have before for stable colors)
    # TODO: remove this if not also used for plot_anno_scores
    if (go_roots_only){
        val_roots = unique(res[[1]][res[[1]][,2] %in% go_ids,1])
    } else {
        val_roots = unique(res[[1]][,1])
    }
    root_names_id = term[match(root_names, term[,2]) ,]
    root_names_id = root_names_id[root_names_id[,2] %in% val_roots, ]
    root_ids = root_names_id[,4]
    
    # add root nodes
    if (test != "contingency"){
        go_ids = c(go_ids, root_ids)
    }
    
    # add root_ids to GO-ids for combining them with output
    matched_root_name = get_names(go_ids, term) # get names also returns root-node-name
    matched_root_name$root_id = root_names_id[match(matched_root_name[,3], root_names_id[,2]), 4]
    
    # check if background genes are defined (optional for hyper)
    if (test == "hyper" & all(in_genes[,2] == 1)){
        genes = NULL
    } else {
        genes = in_genes[,1]
    }
    ### get annotation for nodes
    # TODO: move this check to plotting / refinement function separately - or better: get_anno_genes 
    anno_db = res[[3]][1,2]
    if (anno_db == "custom" && is.null(annotations)){
        stop("Apparently go_enrich was run with custom annotations. Please provide those custom annotations to plot_anno_scores, too. See ?plot_anno_scores")
    }
    
    message("\nGet annotations for categories:")
    anno = get_anno_genes(go_ids, database=anno_db, genes, annotations, term, graph_path)
    if (is.null(anno)) return(invisible(anno)) # no annotations - warning from get_anno_genes
    # add root_ids and scores to nodes
    anno = cbind(anno, root_id=matched_root_name[match(anno[,1], matched_root_name[,1]), "root_id"])
    anno$root_id = as.character(anno$root_id)
    anno_scores = cbind(anno, in_genes[match(anno[,2], in_genes[,1]), 2:ncol(in_genes)])
    colnames(anno_scores)[4:ncol(anno_scores)] = colnames(in_genes)[2:ncol(in_genes)]
    
    # for default background, bg-genes are not in res[[2]], match yields NA, convert to 0
    if (test == "hyper"){
        anno_scores[is.na(anno_scores[,4]), 4] = 0 # default 0 for hyper (NA if background not defined)
    }
    
    return(anno_scores)
}
    
    
