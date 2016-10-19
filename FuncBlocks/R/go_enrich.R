
# run original "FUNC" in subdirectory ("tmp")
# reads term.txt, term2term.txt, graph_path.txt from /r1/steffi_grote/R_packages/term_tables 

# argument "blocks": if True, read gene_coordinates, modify "root-node", and call hyper_randset_blocks

### zum Testen:
#setwd("~/ownCloud/forAkeyPaper/FUNCpackage_Akey")
#test_genes = paste(rep('FABP', 5), c(2,4:7), sep='')
#bg_genes = c('NCAPG', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1')
#genes = c(rep(1,length(test_genes)), rep(0,length(bg_genes)))
#names(genes) = c(test_genes, bg_genes)
#test="hyper"
#blocks=FALSE
#n_randsets=10

go_enrich=function(genes, test="hyper", blocks=FALSE, n_randsets=1000)
{

	if (blocks==TRUE){
		print(head(gene_coords_bed))
		stop("This option is currently not available :(. But gene-positions are already included! :)") 
#		gene_coords = read.table("~/ownCloud/forAkeyPaper/gene_coordinates.bed")[,1:4]
#		print(head(gene_coords))
#		if (!(all(names(genes) %in% gene_coords[,4]))){ 
#			stop("Not all genes have coordinates in 'gene_coordinates.bed'!")}
	}
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
		n_randsets=round(n_randsets)
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
	} else (stop("Not a valid test. Please use 'hyper' or 'wilcoxon'."))

	# Create input and output folder and write file with missing genes
	directory = tempdir()
#	dir.create("./tmp"); directory = "./tmp"
	
	message("get GOs for genes")
#	go = read.table("/mnt/expressions/miguel/neandertal_int/mart_export.txt",as.is=T,sep="\t",head=T)
	go = go_anno[go_anno[,2] %in% names(genes) & go_anno[,1]!="",2:1]  
	go$value = genes[match(go[,1], names(genes))]
	# subset to GOs present in term.txt (which is now included in sys.data)
#	term = read.table("/r1/people/steffi_grote/R_packages/term_tables/term.txt" ,sep="\t", quote="", comment.char="", as.is=TRUE)
	# remove obsolete terms
	term = term[term[,5]==0,]
	go = go[go[,2] %in% term[,4],]
	
	# NEW: write ontolgy-graph tables to tmp-directory (which are now included in sys.data)
	write.table(term,file=paste(directory, "/term.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	write.table(term2term,file=paste(directory, "/term2term.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	write.table(graph_path,file=paste(directory, "/graph_path.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
	
	## loop over different root-nodes
	# nur ein "universal root" hat die root-flag in term.txt
	root_node_ids= c("GO:0003674","GO:0008150","GO:0005575")
	root_nodes = term[match(root_node_ids, term[,4]),2]
	out = data.frame()

	for (r in 1:length(root_nodes)){
		root_node = root_nodes[r]
		message(root_node)
		# subset the input data to GOs that belong to current root node (col 3 in term.txt is root-node)
		root_data = go[term[match(go[,2],term[,4]),3]==root_node,]
		
		# prepare root_data data (infile-data and root like in separate_taxonomies.pl)
		# paste annotations
		xx = tapply(root_data[,2],root_data[,1],function(x){paste(x,collapse=" ")})
		
		# create "root" dataframe 
		if(blocks){
			# Fuer ABAEnrichmentBlocks: add gene-coords to "root-node" file
			gene = as.character(names(xx))
			gene_position = gene_coords[match(gene, gene_coords[,4]),1:3]	
			root = data.frame(genes=gene, gene_position ,goterms=as.character(xx))
		} else {
			root=data.frame(genes=as.character(names(xx)),goterms=as.character(xx)) 
		}

		# "infile-data" (wilcox) / "Allen:4005-changed" (hyper)
		if (test=="hyper"){
			# subset to test genes
			root_data=root_data[root_data[,3]==1,]	
			yy=tapply(root_data[,2],root_data[,1],function(x){paste(x,collapse=" ")})
			# create infile-data dataframe
			infile_data=data.frame(genes=as.character(names(yy)),goterms=as.character(yy))
		} else if(test=="wilcoxon"){
			infile_data = unique(root_data[,c(1,3)])
		}

		# write root_data-files to tmp-directory
		write.table(infile_data,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste(directory,"/infile-data",sep=""))
		write.table(root,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste(directory,"/",root_node_ids[r],sep=""))
		
		# for debugging: save input-files
#		system(paste("cp ",directory,"/infile-data ", directory, "/infile_data_", root_node, sep=""))
							
		##	3b) run func
#			stop("Don't run FUNC!")
		if(test=="hyper"){
			if (blocks){
				run_func(hyper_randset_blocks, hyper_category_test, directory, root_node_ids[r], n_randsets)	
			} else {
				run_func(hyper_randset, hyper_category_test, directory, root_node_ids[r], n_randsets)	
			}		
		} else if (test=="wilcoxon"){				
			run_func(wilcox_randset, wilcox_category_test, directory, root_node_ids[r], n_randsets)
		}	
	
		groupy=read.table(paste(directory,"/category_test_out",sep=""))
		system(paste("mv ", directory, "/category_test_out ",directory, "/out_", root_node, sep=""))
		out = rbind(out, groupy)
	}
	
	# add GO-names and sort
	namen = term[match(out[,1],term[,4]),2:3]
	out = data.frame(namen[,2],out[,1], namen[,1], out[,2:ncol(out)])
	out = out[order(out[,7], out[,5]),]
	
	if (test == "hyper"){
		colnames(out)=c("ontology","node_id","node_name","raw_p_underrep","raw_p_overrep","FWER_underrep","FWER_overrep")
	} else if (test == "wilcoxon"){
		colnames(out)=c("ontology","node_id","node_name","raw_p_low_rank","raw_p_high_rank","FWER_low_rank","FWER_high_rank")
	}
	
	message("\nDone.")
	return(out)	
}	

