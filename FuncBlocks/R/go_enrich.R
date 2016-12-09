# run "FUNC" with default GO and GO-annotations (update once in a while...) TODO: also allow user input
# wilcoxon rank test or hypergeometric test

### Arguments
# genes: vector of 0/1 (hypergeometric) or float (wilcoxon) with gene identifiers as names (HGNC-symbol), or chromosomal regions chr:start-stop (only for hypergeometric)
# test: "hyper" or "wilcoxon"
# n_randsets: number of random-sets
# gene_len: randomset is dependent on length of genes
# circ_chrom: for regions input: random regions are on same chrom and allowed to overlap multiple bg-regions

### Testing:
#/r1/people/steffi_grote/R_packages/FuncBlocks_package/test.R)

go_enrich=function(genes, test="hyper", n_randsets=1000, gene_len=FALSE, circ_chrom=FALSE)
{
	
	ref_genome = "grch37" # TODO: add grch38 as parameter and gene_coords_grch38 to sysdata.R 
	
	#####	1. Check arguments and define parameters
	
	## Check arguments
	# general
	message("Checking arguments...")
	if (length(genes)==0){
		stop("Please enter genes.")
	}	
	if (length(names(genes))==0){
		stop("Please add gene identifiers as names to 'genes' vector.")
	}
	if(length(n_randsets)>1 || !is.numeric(n_randsets) || n_randsets<1){
		stop("Please define 'n_randsets' as a positive integer.")
	}
	if(n_randsets != round(n_randsets)){
		n_randsets = round(n_randsets)
		warning(paste("'n_randsets' is expected to be an integer and was rounded to ",n_randsets,".",sep=""))
	}
	# test-specific arguments
	if (test=="hyper"){
		if(!all(genes %in% c(0,1))){
			stop("Not a valid 'genes' argument for hypergeometric test. Please use a vector of 0/1.")	
		}	
		if(sum(genes)==0){
			stop("Only 0 in genes vector. Please enter test genes.")	
		}
	} else	if (test=="wilcoxon"){
		if(!is.numeric(genes)){
			stop("Not a valid 'genes' argument. Please use a numeric vector.")	
		}	
		if(gene_len == TRUE){
			stop("Argument 'gene_len = TRUE' can only be used with 'test = 'hyper''.")
		}	
	} else (stop("Not a valid test. Please use 'hyper' or 'wilcoxon'."))
	
	
	#####	2. Prepare for FUNC	
	
	# Create input and output folder and write file with missing genes
	directory = tempdir()
#	dir.create("./tmp"); directory = "./tmp"

	# load gene coordinates
	gene_coords = get(paste("gene_coords_", ref_genome, sep=""))
	
	# are genes or genomic regions ('blocks') given as input?
	blocks = FALSE
	identifier = detect_identifier(names(genes)[1]) # "blocks" or "hgnc_symbol"
	if (identifier=="blocks"){
#		identifier = "hgnc_symbol" # gene-name 
		blocks = TRUE
		# check that background region is specified
		if(sum(genes) == length(genes)){
			stop("All values of the 'genes' input are 1. Using chromosomal regions as input requires defining background regions with 0.")
		}
		# check that test is hyper
		if (test != "hyper"){
			stop("chromosomal regions can only be used with test='hyper'.")
		}	
		# warn if gene_len=TRUE, although regions are used
		if (gene_len == TRUE){
			warning("Unused argument: 'gene_len = TRUE'.")
		}		
		# convert coords from genes-vector to bed format, SORT, CHECK and extract genes overlapping regions
		regions = get_genes_from_regions(genes, gene_coords, circ_chrom) # gene_coords from sysdata.rda HGNC
		test_regions = regions[[1]]		
		bg_regions = regions[[2]]
		genes = regions[[3]]
		
		# avoid scientific notation in regions (read in c++)
		test_regions = format(test_regions, scientific=FALSE, trim=TRUE)
		bg_regions = format(bg_regions, scientific=FALSE, trim=TRUE)

		message("Candidate regions:")
		print(test_regions)
		message("Background regions:")
		print(bg_regions)		
	
		# write regions to files
		write.table(test_regions,file=paste(directory, "/test_regions.bed",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
		write.table(bg_regions,file=paste(directory, "/bg_regions.bed",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")		
	} else if (circ_chrom == TRUE){
		# warn if circ_chrom=TRUE, although individual genes are used
		warning("Unused argument: 'circ_chrom = TRUE'.")
	} 	
	
	message("get GOs for genes...")
	# remove obsolete terms
	term = term[term[,5]==0,]
	# subset to GOs present in term.txt and flip colums
	go = go_anno[go_anno[,1] %in% term[,4] & go_anno[,1]!="", 2:1]
	
	# check which genes dont have GOs annotated (GOs contained in the ontology)
	remaining_genes = genes[names(genes) %in% go[,1]] # restrict
	not_in = unique(names(genes)[!(names(genes) %in% names(remaining_genes))]) # removed	
	if (length(not_in) > 0 && !blocks){	# this message is usually too long when blocks are used. 
		not_in_string = paste(not_in,collapse=", ")
		warning(paste("No GO-annotation for genes: ",not_in_string,".\n  These genes were not included in the analysis.",sep=""))
	}
	# restrict to genes that have coordinates and warn about the rest
	if (gene_len){
		no_coords = unique(names(remaining_genes)[!(names(remaining_genes) %in% gene_coords[,4])]) # removed
		if (length(no_coords) > 0){
			remaining_genes = remaining_genes[!(names(remaining_genes) %in% no_coords)] # restrict
			no_coords_string = paste(no_coords,collapse=", ")
			warning(paste("No coordinates available for genes: ",no_coords_string,".\n  These genes were not included in the analysis.",sep=""))
			not_in = c(not_in, no_coords)
		}
	}
	# after removing genes without expression data or coordinates: are enough genes left?
	if (length(not_in) > 0){
		# is any gene in the data? if not, stop.
		if (length(remaining_genes) == 0) {
			stop("No requested genes in data.")
		}
		# for 0/1 data: are test genes and background genes in the data (background only if background specified)?
		if (test=="hyper" & sum(remaining_genes)==0) {
			stop("No requested test genes in data.")
		}
		if (test=="hyper" & sum(genes)!=length(genes) & sum(remaining_genes)==length(remaining_genes)) {
			stop("No requested background genes in data.")
		}
		if (test=="wilcoxon" & length(remaining_genes) < 2) {
			stop(paste("Less than 2 genes have annotated GOs.",sep=""))
		}
	}
			
	
	# subset to input genes (unless test=hyper & no background genes defined)
	if(!(test=="hyper" & sum(genes) == length(genes))){
		go = go[go[,1] %in% names(remaining_genes),] 
	}
	# add value for genes (1/0 for hyper, scores for wilcox) 
	go$value = genes[match(go[,1], names(genes))]	


	# write ontolgy-graph tables to tmp-directory (included in sysdata.rda)
	write.table(term,file=paste(directory, "/term.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	write.table(term2term,file=paste(directory, "/term2term.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	write.table(graph_path,file=paste(directory, "/graph_path.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")	


	#####	3. Loop over GO root-nodes and run FUNC
	
	## loop over different root-nodes
	# only "universal root" has root-flag in term.txt -> define explicitly
	root_node_ids= c("GO:0003674","GO:0008150","GO:0005575")
	root_nodes = term[match(root_node_ids, term[,4]),2]
	out = data.frame()

	for (r in 1:length(root_nodes)){		
		root_node = root_nodes[r]
		root_id = root_node_ids[r]
		message(paste("\n\nProcessing root node: ", root_node,"...\n", sep=""))
		
		# subset the input data to GOs that belong to current root node (col 3 in term.txt is root-node)
		input = go[term[match(go[,2],term[,4]),3]==root_node,]
		
		# prepare input data ('infile-data' and 'root' like in separate_taxonomies.pl in FUNC)
		# "infile-data": one column with test genes (hyper) OR genes and associated scores (wilcox) 
		if (test=="hyper"){
			# subset to test genes
			infile_data = data.frame(genes = unique(input[input[,3] == 1,1]))
		} else if(test=="wilcoxon"){
			infile_data = unique(input[,c(1,3)])
		}	
			
		# create "root" dataframe (all genes (test and bg) with GO-annotations, file named with root_node_id)
		# add gene-coordinates, despite for classic FUNC option 
		if (blocks || gene_len){
			# one line per gene: gene | chrom | start | end | GO1 GO2 GO3
			go_string=tapply(input[,2],input[,1],function(x){paste(x,collapse=" ")}) # paste annotations
			gene = as.character(names(go_string))
			# add coordinates
			gene_position = gene_coords[match(gene, gene_coords[,4]),1:3]	
			root = data.frame(genes=gene, gene_position ,goterms=as.character(go_string))		
			# remove genes with unknown coordinates (possible for default bg in gene_len-option) 
			if (gene_len){					
				root = root[!is.na(root[,3]),] 			
			}	
			# NEW just in case - avoid scientific notation of gene coordinates	
			root[,3:4] = format(root[,3:4], scientific=FALSE, trim=TRUE)		
		} else {
			# one line per gene: gene | GO1 GO2 GO3
			go_string = tapply(input[,2],input[,1],function(x){paste(x,collapse=" ")}) # paste annotations
			root = data.frame(genes=as.character(names(go_string)),goterms=as.character(go_string))
		}
	
		# write root_data-files to tmp-directory
		write.table(infile_data,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste(directory,"/infile-data",sep=""))
		write.table(root,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste(directory,"/",root_id,sep=""))
		
		# for debugging: save input-files
#		system(paste("cp ",directory,"/infile-data ", directory, "/infile_data_", root_id, sep=""))
		
		message("Run Func...\n")					
		if (test=="hyper"){
			# randomset
			if (blocks & circ_chrom){
				hyper_randset(paste(directory,"/",root_id,sep=""), n_randsets, directory, root_id, "roll")
			} else if (blocks){
				hyper_randset(paste(directory,"/",root_id,sep=""), n_randsets, directory, root_id,"block")				
			} else if (gene_len){						
				hyper_randset(paste(directory,"/",root_id,sep=""), n_randsets, directory, root_id, "gene_len")	
			} else {						
				hyper_randset(paste(directory,"/",root_id,sep=""), n_randsets, directory, root_id, "classic")
			}	
#			stop("only randomset")
			# category test
			hyper_category_test(paste(directory, "/randset_out",sep=""), paste(directory,"/category_test_out", sep=""), 1, root_id)
		} else if (test=="wilcoxon"){
			wilcox_randset(paste(directory,"/",root_id,sep=""), n_randsets, directory, root_id)
			# category test
			wilcox_category_test(paste(directory, "/randset_out",sep=""), paste(directory,"/category_test_out", sep=""), 1, root_id)
		}

		# read Func output
		groupy=read.table(paste(directory,"/category_test_out",sep=""))
		# for debugging: save output-files per root-node
		#system(paste("mv ", directory, "/category_test_out ",directory, "/out_", root_id, sep=""))
			
		# NEW check that FWER order follows p-value order (per root_node)
		colnames(groupy)=c("node_id","p_under","p_over","FWER_under","FWER_over")
		groupy_sorted = groupy[signif(round(groupy$p_over,12), -groupy$FWER_over),]
		if(any(groupy_sorted$FWER_over != cummax(groupy_sorted$FWER_over))){
			print(data.frame(groupy_sorted[,c(1,3,5)], FWER_check=groupy_sorted$FWER_over == cummax(groupy_sorted$FWER_over)))
			stop("FWER_over does not strictly follow p_over. This looks like a bug.\n  Please contact steffi_grote@eva.mpg.de or david_reher@eva.mpg.de.")
		}
		groupy_sorted = groupy[order(signif(groupy$p_under,12), -groupy$FWER_under),]	
		if(any(groupy_sorted$FWER_under != cummax(groupy_sorted$FWER_under))){
			print(data.frame(groupy_sorted[,c(1,2,4)], FWER_check=groupy_sorted$FWER_under == cummax(groupy_sorted$FWER_under)))
			stop("FWER_under does not strictly follow p_under. This looks like a bug.\n  Please contact steffi_grote@eva.mpg.de or david_reher@eva.mpg.de.")
		}
				
		out = rbind(out, groupy)
	} # end root_nodes
	
	# add GO-names and sort
	namen = term[match(out[,1],term[,4]),2:3]
	out = data.frame(namen[,2],out[,1], namen[,1], out[,2:ncol(out)])
	out = out[order(out[,7], out[,5]),]
	
	if (test == "hyper"){
		colnames(out)=c("ontology","node_id","node_name","raw_p_underrep","raw_p_overrep","FWER_underrep","FWER_overrep", "n_candidate_expected", "n_candidate_real")
	} else if (test == "wilcoxon"){
		colnames(out)=c("ontology","node_id","node_name","raw_p_low_rank","raw_p_high_rank","FWER_low_rank","FWER_high_rank")
	}
	# also return input genes (reduced to those with expression data, candidate genes(no bg defined), with coords(gene_len==T))
	final_output = list(results=out, genes=remaining_genes)
	message("\nDone.")

	return(final_output)
}	

