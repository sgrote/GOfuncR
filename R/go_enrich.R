
# run "FUNC" with integrated GO-graph and GO-annotations from OrganismDb or OrgDb packages
# wilcoxon rank test, hypergeometric test, binomial test, 2x2-contingency-test

go_enrich=function(genes, test="hyper", n_randsets=1000, organismDb="Homo.sapiens", gene_len=FALSE, regions=FALSE, circ_chrom=FALSE, silent=FALSE, domains=NULL, orgDb=NULL, txDb=NULL, annotations=NULL, gene_coords=NULL, godir=NULL)
{
    
    #####   1. Check arguments and define parameters
    
    ## still allow vector as input for hyper and wilcox (like in older versions)
    if (!((is.vector(genes) && !is.null(names(genes))) || is.data.frame(genes))){
        stop("Please provide a data frame as 'genes' input (also named vector is still accepted for hypergeometric or wilcoxon rank-sum test).")
    }
    if (is.vector(genes)){
        genes = data.frame(gene=names(genes), score=unname(genes), stringsAsFactors=FALSE)
    } else {
        genes[,1] = as.character(genes[,1])
    }

    ## Check arguments
    if (!silent) message("Checking arguments...")

    # genes
    if (anyNA(genes)){
        stop("NAs found in 'genes'-input. Please provide finite values only.")
    }
    if (test=="hyper"){
        if (!(is.data.frame(genes) && ncol(genes)==2 && is.numeric(genes[,2]))){
            stop("Please provide a data frame with 2 columns [gene, 0/1] as input for hypergeometric test.")
        }
        if (!(all(genes[,2] %in% c(1,0)))){
            stop("Please provide only 1/0-values in 2nd column of 'genes'-input for hypergeometric test.")
        }
        if (sum(genes[,2])==0){
            stop("Only background genes defined (only 0s in 'genes[,2]'). Please provide candidate genes denoted by 1.")
        }
    } else  if (test=="wilcoxon"){
        if (!(is.data.frame(genes) && ncol(genes)==2  && is.numeric(genes[,2]))){
            stop("Please provide a data frame with 2 columns [gene, score] as input for wilcoxon rank-sum test.")
        }
        if (nrow(genes) < 2){
            stop("Only one gene provided as input.")
        }
    } else if (test=="binomial"){
        if (!(is.data.frame(genes) && ncol(genes)==3 && all(sapply(genes,is.numeric) == c(0,1,1)))){
            stop("Please provide a data frame with columns [gene, count1, count2] as input for binomial test.")
        }
        if (any(genes[,2:3] < 0) | any(genes[,2:3] != round(genes[,2:3]))){
            stop("Please provide non-negative integers in columns 2-3.")
        }
    } else if (test=="contingency"){
        if (!(is.data.frame(genes) && ncol(genes)==5 && all(sapply(genes,is.numeric) == c(0,1,1,1,1)))){
            stop("Please provide a data frame with columns [gene, count1A, count2A, count1B, count2B] as input for contingency table test.")
        }
        if (any(genes[,2:5] < 0) | any(genes[,2:5] != round(genes[,2:5]))){
            stop("Please provide non-negative integers in columns 2-5.")
        }
    } else {
        stop("Not a valid test. Please use 'hyper', 'wilcoxon', 'binomial' or 'contingency'.")
    }

    # check for multiple assignment for one gene
    genes = unique(genes) # allow for multiple assignment of same value
    multi_genes = sort(unique(genes[,1][duplicated(genes[,1])]))
    if (length(multi_genes) > 0){
        if (test == "hyper"){
            # allow but remove background definition of candidate genes
            # (FUNC just takes one value per gene, candidate genes are implictly part of background)
            candi_genes = genes[genes[,2]==1, 1]
            genes = genes[!(genes[,1] %in% candi_genes & genes[,2] == 0),]
        } else {
            stop("Genes with multiple assignment in input: ", paste(multi_genes,collapse=", "))
        }
    }

    # other arguments
    if (length(n_randsets)>1 || !is.numeric(n_randsets) || n_randsets<1){
        stop("Please define 'n_randsets' as a positive integer.")
    }
    if (n_randsets != round(n_randsets)){
        n_randsets = round(n_randsets)
        warning("'n_randsets' is expected to be an integer and was rounded to ",n_randsets,".")
    }
    if (!is.logical(gene_len)){
        stop("Please set gene_len to TRUE or FALSE.")
    }
    if (!is.logical(circ_chrom)){
        stop("Please set circ_chrom to TRUE or FALSE.")
    }
    if (gene_len & !(test == "hyper")){
        stop("Argument 'gene_len = TRUE' can only be used with 'test = 'hyper''.")
    }
    if (regions & !(test == "hyper")){
        stop("Argument 'regions = TRUE' can only be used with 'test = 'hyper''.")
    }
    if (circ_chrom & !regions){
        stop("Argument 'circ_chrom = TRUE' can only be used with 'regions = TRUE'.")
    }
    if (!is.null(gene_coords)){
        if (test != "hyper"){
            stop("Argument 'gene_coords' can only be used with 'test = 'hyper''.")
        }
        if (ncol(gene_coords) != 4){
            stop("'gene_coords' has to be a data.frame with 4 columns: [gene, chromosome, start, end].")
        }
        # reorder columns like for db-gene-coords
        gene_coords = gene_coords[,c(2:4,1)]
        # add 'chr' if missing (coord_db has it too)
        if (!startsWith(as.character(gene_coords[1,1]), "chr")){ 
            gene_coords[,1] = paste0("chr", gene_coords[,1])
        }
        # check for multiple assignment of coords to genes
        gene_coords = unique(gene_coords)
        multi_coords = sort(unique(gene_coords[,4][duplicated(gene_coords[,4])]))
        if (length(multi_coords) > 0){
            stop("Genes with multiple coordinates in 'gene_coords': ", paste(multi_coords,collapse=", "))
        }
    }
    
    # evaluate database + go + opt combinations and get table including versions
    databases = eval_db_input(organismDb, godir, orgDb, annotations, txDb, regions, gene_len, gene_coords)
    anno_db = databases[databases[,1] == "go_annotations", 2]
    coord_db = databases[databases[,1] == "gene_coordinates", 2] # might be NA
    entrez_db = databases[databases[,1] == "symbol_to_entrez", 2] # might be NA
    # remove NA-entries for output
    databases = databases[!(is.na(databases[,2])),]
    row.names(databases) = 1:nrow(databases)
    
    
    #####   2. Prepare for FUNC
    
    # GO-graph
    # Create tempfile prefix (in contrast to tempdir() alone, this allows parallel processing)
    directory = tempfile()
    # unless custom GO defined, write ontolgy-graph tables to tmp-directory (included in sysdata.rda)
    if (is.null(godir)){
        if (!silent) message("Write temporary files...")
        go_paths = paste0(directory, c("_term.txt", "_term2term.txt", "_graph_path.txt"))
        write.table(term, go_paths[1], col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")
        write.table(term2term, go_paths[2], col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
        write.table(graph_path, go_paths[3], col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
    } else {
        # ~ seems to not work within C++ code
        if (substring(godir, 1, 1) == "~"){
            home = Sys.getenv("HOME")
            godir = paste0(home, substring(godir, 2))
        }
        go_paths = paste0(godir, c("/term.txt", "/term2term.txt", "/graph_path.txt"))
        term = read.table(go_paths[1], sep="\t", quote="", comment.char="", as.is=TRUE)
    }
    
    # convert genomic regions to single genes (also check them and write to file for C++)
    if (regions){
        if (!silent) message("get genes from regions...")
        block_info = blocks_to_genes(directory, genes, coord_db, entrez_db, gene_coords, circ_chrom, silent)
        genes = block_info[[1]] # 'classic' test and bg genes
        gene_coords = block_info[[2]] # gene coords for test and bg genes
    }
    
    ### get GO-annotations
    # if test=hyper and default background get annotations for all genes in database
    if (test=="hyper" && all(genes[,2]==1)){
        go_anno = get_anno_categories(database=anno_db, annotations=annotations, term_df=term, silent=silent)
    } else {
        go_anno = get_anno_categories(genes[,1], database=anno_db, annotations=annotations, term_df=term, silent=silent)
    }
    
    ### get root-nodes from annotations + term 
    # (custom ontology with different root_nodes does not need to be defined in 'domains')
    root_nodes = sort(unique(term[term[,4] %in% go_anno[,2],3]))
    # remove empty root node (very few cases)
    root_nodes = root_nodes[root_nodes != ""]
    if (!silent) message("Found root_nodes: ", paste(root_nodes, collapse=", "))
    if (!is.null(domains)){
        if (!(all(domains %in% root_nodes))){
            stop("'domains' must be in the set of '", paste(root_nodes, collapse=", "),"'.")
        }
        root_nodes = domains
    }    
    # check that all root nodes in term
    root_node_ids = term[match(root_nodes, term[,2]),4]
    if (!(is.null(godir)) && any(is.na(root_node_ids))){
        stop("Not all root_nodes present in term.txt.")
    }
    
    # subset genes to annotated genes (also reduced to internal ontology)
    if (!silent) message("Remove invalid genes...")
    gene_values = genes[genes[,1] %in% go_anno[,1],] # restrict
    if (nrow(gene_values) == 0) {
        stop("No GO-annotations for input genes.") 
    }
    not_in = genes[,1][!genes[,1] %in% gene_values[,1]] # removed
    if (length(not_in) > 0 && !regions){ # this message is usually too long when regions are used
        not_in_string = paste(not_in,collapse=", ")
        warning("No GO-annotation for genes: ",not_in_string,".\n  These genes were not included in the analysis.")
    }

    ### get coordinates
    if (gene_len & is.null(gene_coords)){
        # load gene coordinates
        if (test=="hyper" && all(genes[,2]==1)){
            gene_coords = suppressWarnings(get_all_coords(coord_db, entrez_db, silent))
        } else {
            gene_coords = suppressWarnings(get_gene_coords(genes[,1], coord_db, entrez_db, silent))
        }
    }
    # subset to genes that have coordinates and warn about the rest
    if (gene_len){ # only for test==hyper
        no_coords = gene_values[!(gene_values[,1] %in% gene_coords[,4]), 1] # removed
        if (length(no_coords) > 0){
            gene_values = gene_values[!(gene_values[,1] %in% no_coords),] # restrict
            no_coords_string = paste(no_coords,collapse=", ")
            warning("No coordinates available for genes: ",no_coords_string,".\n  These genes were not included in the analysis.")
            not_in = c(not_in, no_coords)
        }
    }
    
    # after removing genes without annotations or coordinates: are enough genes left?
    if (length(not_in) > 0){
        # is any gene in the data? if not, stop.
        if (nrow(gene_values) == 0) {
            stop("No requested genes in data.")
        }
        # for 0/1 data: are test genes and background genes in the data (background only if background specified)?
        if (test=="hyper" && sum(gene_values[,2])==0) {
            stop("No requested candidate genes in data.")
        }
        if (test=="hyper" && 0 %in% genes[,2] && all(gene_values[,2]==1)) {
            stop("No requested background genes in data.")
        }
        # at least two for wilcoxon
        if (test=="wilcoxon" && nrow(gene_values) < 2) {
            stop("Less than 2 genes have annotated GO-categories.")
        }
    }
    
    
    #####   3. Loop over GO root-nodes and run FUNC
    
    out = list()
    min_p = list()
    
    for (r in seq_along(root_nodes)){
        root_node = root_nodes[r]
        root_id = root_node_ids[r]
        if (!silent) message("\n\nProcessing root node: ", root_node,"...\n")
        
        # subset GO-annotations to nodes that belong to current root node (col 3 in term.txt is root-node)
        gene_go_root = go_anno[term[match(go_anno[,2],term[,4]),3]==root_node,]
        
        # skip root-node if no annotations of input genes (at least two for wilcoxon) 
        if (nrow(gene_go_root) == 0 || (test == "wilcoxon" && nrow(gene_go_root) < 2)){
            warning("No GO-annotations for root node '", root_node,"'.")
            message("\nSkipping root node '", root_node,"'.\n")
            next
        }

        ### prepare input data ('infile-data' and 'root' in c++ scripts)

        # 'root'
        # gene | (chrom | start | end) | GO1 GO2 GO3
        go_string = tapply(gene_go_root[,2],gene_go_root[,1],function(x){paste(x,collapse=" ")})
        gene = as.character(names(go_string))
        # add gene-coordinates, despite for classic FUNC options
        if (regions || gene_len){
            # one line per gene: gene | chrom | start | end | GO1 GO2 GO3
            # add coordinates
            gene_position = gene_coords[match(gene, gene_coords[,4]),1:3]   
            root = data.frame(genes=gene, gene_position ,goterms=as.character(go_string))       
            # remove genes with unknown coordinates (possible for default bg in gene_len-option) 
            if (gene_len){                  
                root = root[!is.na(root[,3]),]          
            }
            # just in case - avoid scientific notation of gene coordinates  
            root[,3:4] = format(root[,3:4], scientific=FALSE, trim=TRUE)        
        } else {
            # one line per gene: gene | GO1 GO2 GO3
            root = data.frame(genes=gene, goterms=as.character(go_string))
        }

        # 'infile-data'
        # gene | value1 | value2 ...
        # subset to genes annotated to current root-node
        infile_data = gene_values[gene_values[,1] %in% root[,1],]
        # for hyper infile-data is just a column with canidate genes
        if (test=="hyper"){
            infile_data = infile_data[infile_data[,2]==1, 1]
        }
        # write input-files to tmp-directory
        write.table(infile_data,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste(directory,"_infile-data",sep=""))
        write.table(root,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste(directory,"_",root_id,sep=""))
        
        # for debugging: save input-files
#       system(paste("cp ",directory,"_infile-data ", directory, "_infile_data_", root_id, sep=""))
        
        if (!silent) message("Run Func...\n")
        if (test == "hyper"){
            if (regions & circ_chrom){
                hyper_randset(paste0(directory,"_",root_id), n_randsets, directory, go_paths[1], go_paths[2], go_paths[3], root_id, "roll" , silent)
            } else if (regions){
                hyper_randset(paste0(directory,"_",root_id), n_randsets, directory, go_paths[1], go_paths[2], go_paths[3], root_id, "block", silent)
            } else if (gene_len){
                hyper_randset(paste0(directory,"_",root_id), n_randsets, directory, go_paths[1], go_paths[2], go_paths[3], root_id, "gene_len", silent)
            } else {
                hyper_randset(paste0(directory,"_",root_id), n_randsets, directory, go_paths[1], go_paths[2], go_paths[3], root_id, "classic", silent)
            }
            hyper_category_test(directory, 1, root_id, silent)
        } else if (test == "wilcoxon"){
            wilcox_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, go_paths[1], go_paths[2], go_paths[3], root_id, silent)
            wilcox_category_test(directory, 1, root_id, silent)
        } else if (test == "binomial"){
            binom_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, go_paths[1], go_paths[2], go_paths[3], root_id, silent)
            binom_category_test(directory, 1, root_id, silent)          
        } else if (test == "contingency"){
            conti_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, go_paths[1], go_paths[2], go_paths[3], root_id, silent)
            conti_category_test(directory, 1, root_id, silent)
        }
        
        # read Func output
        groupy = read.table(paste(directory,"_category_test_out",sep=""))
        # for debugging: save output-files per root-node
        #system(paste("mv ", directory, "_category_test_out ",directory, "_out_", root_id, sep=""))
        out[[r]] = groupy
        
        # NEW: save minimum-p-values from randomsets (used for refinement)
        min_p_domain = read.table(paste(directory, "_min_p", sep=""))
        min_p[[root_node]] = min_p_domain
        
        
    } # end root_nodes

    out = do.call(rbind, out)
    
    # add GO-names and sort
    namen = term[match(out[,1],term[,4]),2:3]
    out = data.frame(namen[,2],out[,1], namen[,1], out[,2:ncol(out)])
    out[,1:3] = apply(out[,1:3], 2, as.character)
    # switch columns for binomial test high B - then high A to be more consistent
    if (test == "binomial"){
        out[,4:7] = out[,c(5,4,7,6)]
    }
    # also sort on FWER_underrep and raw_p_underrep and GO-ID (ties in tail) 
    out = out[order(out[,7], out[,5], -out[,6], -out[,4], out[,1], out[,2]),] 
    rownames(out) = 1:nrow(out)

    if (test == "hyper"){
        out = out[,1:7]
        colnames(out)=c("ontology","node_id","node_name","raw_p_underrep","raw_p_overrep","FWER_underrep","FWER_overrep")
    } else if (test == "wilcoxon"){
        out = out[,1:7]
        colnames(out)=c("ontology","node_id","node_name","raw_p_low_rank","raw_p_high_rank","FWER_low_rank","FWER_high_rank")
    } else if (test == "binomial"){
        colnames(out)=c("ontology","node_id","node_name","raw_p_high_B","raw_p_high_A","FWER_high_B","FWER_high_A")
    } else if (test == "contingency"){
        colnames(out)=c("ontology","node_id","node_name","raw_p_high_CD","raw_p_high_AB","FWER_high_CD","FWER_high_AB")
    }
    # also return input genes (reduced to those with annotations, candidate genes(no bg defined), with coords(gene_len/regions))
    gene_values = gene_values[mixedorder(gene_values[,1]),]
    gene_values[,1] = as.character(gene_values[,1])
    rownames(gene_values) = 1:nrow(gene_values)
    
    # also return min_p per domain
    domains = rep(names(min_p), each=n_randsets)
    min_p = do.call(rbind, min_p)
    # switch columns for binomial test high B - then high A to be more consistent
    if (test == "binomial"){
        min_p[,c(1,2)] = min_p[,c(2,1)]
    }
    min_p = data.frame(ontology=domains, lower_tail=min_p[,1], upper_tail=min_p[,2])
    
    final_output = list(results=out, genes=gene_values, databases=databases, min_p=min_p)
    if (!silent) message("\nDone.")
    

    return(final_output)
}

