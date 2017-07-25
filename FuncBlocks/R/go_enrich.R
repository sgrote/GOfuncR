
# run "FUNC" with default GO and GO-annotations (update once in a while...) TODO: also allow user input
# wilcoxon rank test, hypergeometric test or binomial test

go_enrich=function(genes, test="hyper", n_randsets=1000, gene_len=FALSE, circ_chrom=FALSE, ref_genome="grch37", silent=FALSE, domains=NULL)
{
	
	#####	1. Check arguments and define parameters
	
	## still allow vector as input for hyper and wilcox (like in older versions)
	if (!((is.vector(genes) && !is.null(names(genes))) || is.data.frame(genes))){
		stop("Please provide a data frame as 'genes' input (also named vector is still accepted for hypergeometric or wilcoxon rank-sum test).")
	}
	if (is.vector(genes)){
		genes = data.frame(gene=names(genes), score=unname(genes), stringsAsFactors=FALSE)
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
	} else	if (test=="wilcoxon"){
		if (!(is.data.frame(genes) && ncol(genes)==2  && is.numeric(genes[,2]))){
			stop("Please provide a data frame with 2 columns [gene, score] as input for wilcoxon rank-sum test.")
		}
		if (nrow(genes) < 2){
			stop("Only one gene provided as input.")
		}
		if (gene_len == TRUE){
			stop("Argument 'gene_len = TRUE' can only be used with 'test = 'hyper''.")
		}
	} else if (test=="binomial"){
		if (!(is.data.frame(genes) && ncol(genes)==3 && all(sapply(genes,is.numeric) == c(0,1,1)))){
			stop("Please provide a data frame with columns [gene, count1, count2] as input for binomial test.")
		}
		if (any(c(genes[,2],genes[,3]) != round(c(genes[,2],genes[,3]))) || any(c(genes[,2], genes[,3]) < 0)){
			stop("Please provide a data frame with columns [gene, count1, count2] as input for binomial test. count1 and count2 need to be integers >= 0.")
		}
		if (test=="binomial" && all(c(genes[,2],genes[,3]) == 0)) {
			stop(paste("All input values are 0.",sep=""))
		}
		if (gene_len == TRUE){
			stop("Argument 'gene_len = TRUE' can only be used with 'test = 'hyper''.")
		}
	} else if (test=="contingency"){
		if (!(is.data.frame(genes) && ncol(genes)==5 && all(sapply(genes,is.numeric) == c(0,1,1,1,1)))){
			stop("Please provide a data frame with columns [gene, count1A, count2A, count1B, count2B] as input for contingency table test.")
		}
		if (any(sapply(genes[,2:5], function(x) x<0))){
			stop("Please provide non-negative values.")
		}
	} else {
		stop("Not a valid test. Please use 'hyper', 'wilcoxon', 'binomial' or 'contingency'.")
	}

	# check for multiple assignment for one gene
	genes = unique(genes)
	multi_genes = sort(unique(genes[,1][duplicated(genes[,1])]))
	if (length(multi_genes) > 0){
		stop(paste("Genes with multiple assignment in input:", paste(multi_genes,collapse=", ")))
	}

	# other arguments
	if (!(ref_genome %in% c("grch37","grch38","grcm38"))){
		stop("Please use 'ref_genome=grch37', 'ref_genome=grch38' or 'ref_genome=grcm38'")
	}
	if (length(n_randsets)>1 || !is.numeric(n_randsets) || n_randsets<1){
		stop("Please define 'n_randsets' as a positive integer.")
	}
	if (n_randsets != round(n_randsets)){
		n_randsets = round(n_randsets)
		warning(paste("'n_randsets' is expected to be an integer and was rounded to ",n_randsets,".",sep=""))
	}
	if (!is.logical(gene_len)){
		stop("Please set gene_len to TRUE or FALSE.")
	}
	if (!is.logical(circ_chrom)){
		stop("Please set circ_chrom to TRUE or FALSE.")
	}
	root_nodes = c("molecular_function","biological_process","cellular_component")
	if (!is.null(domains)){
		if (!all(domains %in% root_nodes)){
			stop("'domains' must be in the set of '", paste(root_nodes, collapse=", "),"'.")
		}
		root_nodes = domains
	}

	
	#####	2. Prepare for FUNC	
	
	# Create tempfile prefix (in contrast to tempdir() alone, this allows parallel processing)
	directory = tempfile()
#	dir.create("tempdir"); directory = paste("./tempdir/tempfile",Sys.info()["nodename"],sep="_")

	# load gene coordinates
	gene_coords = get(paste("gene_coords_", ref_genome, sep=""))

	# load GO-annotation
	go_anno = get(paste("go_anno_", ref_genome, sep=""))
	
	# detect identifier: are genes or genomic regions ('blocks') given as input?
	blocks = grepl("^[0-9XY]*:[0-9]*-[0-9]*$", genes[1,1]) # TODO: allow more than [0-9XY] as chroms?

	# convert genomic regions to single genes (get_genes_from_regions)
	if (blocks){
		genes = blocks_to_genes(directory, genes, test, gene_len, gene_coords, circ_chrom, silent)
	} else {
		if (circ_chrom == TRUE){
			# warn if circ_chrom=TRUE, although individual genes are used
			warning("Unused argument: 'circ_chrom = TRUE'.")
		}
	}
	
	if (!silent) message("get GOs for genes...")
	# remove obsolete terms
	term = term[term[,5]==0,]
	# subset to GOs present in term.txt and flip colums [gene, GO]
	gene_go = go_anno[go_anno[,1] %in% term[,4] & go_anno[,1]!="", 2:1]

	# subset genes to annotated genes (GOs contained in the ontology)
	if (!silent) message("Remove invalid genes...")
	gene_values = genes[genes[,1] %in% gene_go[,1],] # restrict
	not_in = genes[,1][!genes[,1] %in% gene_values[,1]] # removed
	if (length(not_in) > 0 && !blocks){	# this message is usually too long when blocks are used. 
		not_in_string = paste(not_in,collapse=", ")
		warning(paste("No GO-annotation for genes: ",not_in_string,".\n  These genes were not included in the analysis.",sep=""))
	}

	# subset GO-annotations to input genes (unless test=hyper & no background genes defined)
	if (!(test=="hyper" && all(gene_values[,2]==1))){
		gene_go = gene_go[gene_go[,1] %in% gene_values[,1],]
	}

	# subset to genes that have coordinates and warn about the rest
	if (gene_len){ # only for test==hyper
		no_coords = gene_values[,1][!(gene_values[,1] %in% gene_coords[,4])] # removed
		if (length(no_coords) > 0){
			gene_values = gene_values[!(gene_values[,1] %in% no_coords)] # restrict
			no_coords_string = paste(no_coords,collapse=", ")
			warning(paste("No coordinates available for genes: ",no_coords_string,".\n  These genes were not included in the analysis.",sep=""))
			not_in = c(not_in, no_coords)
		}
	}
	# after removing genes without expression data or coordinates: are enough genes left?
	if (length(not_in) > 0){
		# is any gene in the data? if not, stop.
		if (nrow(gene_values) == 0) {
			stop("No requested genes in data.")
		}
		# for 0/1 data: are test genes and background genes in the data (background only if background specified)?
		if (test=="hyper" && sum(gene_values[,2])==0) {
			stop("No requested test genes in data.")
		}
		if (test=="hyper" && 0 %in% genes[,2] && all(gene_values[,2]==1)) {
			stop("No requested background genes in data.")
		}
		# at least two for wilcoxon
		if (test=="wilcoxon" && nrow(gene_values) < 2) {
			stop(paste("Less than 2 genes have annotated GOs.",sep=""))
		}
	}

	# write ontolgy-graph tables to tmp-directory (included in sysdata.rda)
	if (!silent) message("Write temporary files...")
	write.table(term,file=paste(directory, "_term.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	write.table(term2term,file=paste(directory, "_term2term.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	write.table(graph_path,file=paste(directory, "_graph_path.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")


	#####	3. Loop over GO root-nodes and run FUNC

	root_node_ids = term[match(root_nodes, term[,2]),4]

	## check that all root nodes in term (for custom ontology) 
	# TODO: activate for custom onto
#	if (length(root_node_ids) != length(root_nodes)){
#		stop("Not all root_nodes present in term.txt")
#	}
	
	out = list()

	for (r in 1:length(root_nodes)){		
		root_node = root_nodes[r]
		root_id = root_node_ids[r]
		if (!silent) message(paste("\n\nProcessing root node: ", root_node,"...\n", sep=""))
		
		# subset GO-annotations to nodes that belong to current root node (col 3 in term.txt is root-node)
		gene_go_root = gene_go[term[match(gene_go[,2],term[,4]),3]==root_node,]
		
		# NEW: skip root-node if no annotations of input genes (at least two for wilcoxon) 
		if (nrow(gene_go_root) == 0 || (test == "wilcoxon" && nrow(gene_go_root) < 2)){
			warning(paste0("No GO-annotations for root node '", root_node,"'."))
			message(paste0("\nSkipping root node '", root_node,"'.\n"))
			next
		}

		### prepare input data ('infile-data' and 'root' in c++ scripts)

		# 'root'
		# gene | (chrom | start | end) | GO1 GO2 GO3
		go_string = tapply(gene_go_root[,2],gene_go_root[,1],function(x){paste(x,collapse=" ")})
		gene = as.character(names(go_string))
		# add gene-coordinates, despite for classic FUNC option 
		if (blocks || gene_len){
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
#		system(paste("cp ",directory,"_infile-data ", directory, "_infile_data_", root_id, sep=""))
		
		if (!silent) message("Run Func...\n")
		if (test == "hyper"){
			# randomset
			if (blocks & circ_chrom){
				hyper_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, "roll" , silent)
			} else if (blocks){
				hyper_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id,"block", silent)
			} else if (gene_len){
				hyper_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, "gene_len", silent)
			} else {
				hyper_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, "classic", silent)
			}
			# category test
			hyper_category_test(directory, 1, root_id, silent)
		} else if (test == "wilcoxon"){
			wilcox_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, silent)
			wilcox_category_test(directory, 1, root_id, silent)
		} else if (test == "binomial"){
			binom_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, silent)
			binom_category_test(directory, 1, root_id, silent)			
		} else if (test == "contingency"){
			conti_randset(paste(directory,"_",root_id,sep=""), n_randsets, directory, root_id, silent)
			conti_category_test(directory, 1, root_id, silent)
		}
		
		# read Func output
		groupy = read.table(paste(directory,"_category_test_out",sep=""))
		# for debugging: save output-files per root-node
		#system(paste("mv ", directory, "_category_test_out ",directory, "_out_", root_id, sep=""))
			
		# check that FWER order follows p-value order (per root_node)
		colnames(groupy)[1:5]=c("node_id","p_under","p_over","FWER_under","FWER_over")
		groupy_sorted = groupy[signif (round(groupy$p_over,12), -groupy$FWER_over),]
		if (any(groupy_sorted$FWER_over != cummax(groupy_sorted$FWER_over))){
			print(data.frame(groupy_sorted[,c(1,3,5)], FWER_check=groupy_sorted$FWER_over == cummax(groupy_sorted$FWER_over)))
			stop("FWER_over does not strictly follow p_over. This looks like a bug.\n  Please contact steffi_grote@eva.mpg.de")
		}
		groupy_sorted = groupy[order(signif (groupy$p_under,12), -groupy$FWER_under),]
		if (any(groupy_sorted$FWER_under != cummax(groupy_sorted$FWER_under))){
			print(data.frame(groupy_sorted[,c(1,2,4)], FWER_check=groupy_sorted$FWER_under == cummax(groupy_sorted$FWER_under)))
			stop("FWER_under does not strictly follow p_under. This looks like a bug.\n  Please contact steffi_grote@eva.mpg.de.")
		}
				
		out[[r]] = groupy
	} # end root_nodes

	out = do.call(rbind, out)
	
	# add GO-names and sort
	namen = term[match(out[,1],term[,4]),2:3]
	out = data.frame(namen[,2],out[,1], namen[,1], out[,2:ncol(out)])
	out[,1:3] = apply(out[,1:3], 2, as.character)
	# NEW: also sort on FWER_underrep and raw_p_underrep and GO-ID (ties in tail) 
	out = out[order(out[,7], out[,5], -out[,6], -out[,4], out[,1], out[,2]),] 
	rownames(out) = 1:nrow(out)

	if (test == "hyper"){
		colnames(out)=c("ontology","node_id","node_name","raw_p_underrep","raw_p_overrep","FWER_underrep","FWER_overrep", "n_candidate_expected", "n_candidate_real")
	} else if (test == "wilcoxon"){
		colnames(out)=c("ontology","node_id","node_name","raw_p_low_rank","raw_p_high_rank","FWER_low_rank","FWER_high_rank","ranksum_expected","ranksum_real")
	} else if (test == "binomial"){
		colnames(out)=c("ontology","node_id","node_name","raw_p_high_A","raw_p_high_B","FWER_high_A","FWER_high_B")
	} else if (test == "contingency"){
		colnames(out)=c("ontology","node_id","node_name","raw_p_high_C/D","raw_p_high_A/B","FWER_high_C/D","FWER_high_A/B")
	}
	# also return input genes (reduced to those with expression data, candidate genes(no bg defined), with coords(gene_len==T))
	# NEW: dataframe
	gene_values = gene_values[mixedorder(gene_values[,1]),]

	final_output = list(results=out, genes=gene_values, ref_genome=ref_genome)
	if (!silent) message("\nDone.")

	return(final_output)
}

