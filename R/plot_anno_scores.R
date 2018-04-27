# take a result from go_enrich() and a vector of GO-IDs to plot annotated scores
# plotting type depends on the performed test in go_enrich which is automatically recognized
# (fwer_threshold is not supported anymore, to not discriminate high_A or high_B in binomial, contingency)


plot_anno_scores = function(res, go_ids, annotations=NULL){
    

    ### check input
    # check that res could be go_enrich-output
    if (!(is.list(res) && all(names(res) == c("results","genes","databases")))){
        stop("Please use an object returned from go_enrich as input (list with 3 elements).")
    }
    # check that all go_ids are in res
    if (!all(go_ids %in% res[[1]][,2])){
        inval = go_ids[!go_ids %in% res[[1]][,2]]
        stop("go_ids not present in res (go_enrich result): ", paste(inval,collapse=", "))
    }
    # infer test
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
    
    # check if background genes are defined (optional for hyper)
    if (test == "hyper" & all(in_genes[,2] == 1)){
        genes = NULL
    } else {
        genes = in_genes[,1]
    }
    ### get annotation for nodes
    anno_db = res[[3]][1,2]
    if (anno_db == "custom" && is.null(annotations)){
        stop("Apparently go_enrich was run with custom annotations. Please provide those custom annotations to plot_anno_scores, too. See ?plot_anno_scores")
    }
    anno = get_anno_genes(go_ids, database=anno_db, genes, annotations, term, graph_path)
    if (is.null(anno)) return(invisible(anno)) # no annotations - warning from get_anno_genes
    # add scores to nodes
    anno_scores = cbind(anno, in_genes[match(anno[,2], in_genes[,1]), 2:ncol(in_genes)])

    # aggregate scores in nodes (wilcox: plot score distribution)
    if (test != "wilcoxon"){
        if (test == "hyper"){ 
            # counts of 1 and 0 genes in a node
            anno_scores[is.na(anno_scores[,3]), 3] = 0 # default 0 for hyper (NA if background not defined)
            anno_scores = tapply(anno_scores[,3], anno_scores[,1], function(x) c(sum(x[]), length(x)-sum(x)))
            anno_scores = data.frame(go_id = names(anno_scores), do.call(rbind, anno_scores))
        } else { 
            # sums of scores in a node (binom + conti)
            anno_scores = aggregate(anno_scores[,3:ncol(anno_scores)], list(go_id=anno_scores[,1]), sum)
        }
    }


    ### get annotation for root nodes   (conti independent of root nodes)
    if (test != "contingency"){
        
        # get annotation for root nodes
        root_anno = get_anno_genes(root_ids, database=res[[3]][1,2], genes, annotations, term, graph_path)
        # add scores to root
        root_anno_scores = cbind(root_anno, in_genes[match(root_anno[,2], in_genes[,1]), 2:ncol(in_genes)])
        # order alphabetically (which is rev(order(IDs)) for default GO)
        root_anno_scores = root_anno_scores[rev(order(root_anno_scores[,1])), ]

        # aggregate scores in root nodes
        if (test != "wilcoxon"){
            if (test == "hyper"){ 
                # counts of 1 and 0 genes in a node
                root_anno_scores[is.na(root_anno_scores[,3]), 3] = 0 
                root_anno_scores = tapply(root_anno_scores[,3], root_anno_scores[,1], function(x) c(sum(x[]), length(x)-sum(x)))
                root_anno_scores = data.frame(go_id = names(root_anno_scores), do.call(rbind, root_anno_scores))
            } else { 
                # sums of scores in a node (binom + conti)
                root_anno_scores = aggregate(root_anno_scores[,3:ncol(root_anno_scores)], list(go_id=root_anno_scores[,1]), sum)
            }
            # add colors and root_node_name
            root_anno_scores$root_name = get_names(root_anno_scores[,1], term)[,2]
            root_anno_scores$root_col = root_cols[match(root_anno_scores[,1], root_cols[,1]), 2]
            # merge nodes with root node info
            matched_root_name = get_names(anno_scores[,1], term)[,3] # get names
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
    
    
    

