## run FUNC with either "hyper" or "wilcoxon"
#(input are files except for n_randsets)
#(output is file)

run_func=function(FUN1, FUN2, directory, root_node, n_randsets)
{
	message("run_func...")
	# randset
	FUN1(
		paste(directory,"/",root_node,sep=""),
		paste(directory,"/infile-data",sep=""),
		n_randsets,
		paste(directory, "/randset_out",sep=""),
		paste(directory, "/term.txt",sep=""),
		paste(directory, "/graph_path.txt",sep=""),
		paste(directory, "/term2term.txt",sep=""),
#		"/r1/people/steffi_grote/R_packages/term_tables/term.txt", 
#		"/r1/people/steffi_grote/R_packages/term_tables/graph_path.txt",
#		"/r1/people/steffi_grote/R_packages/term_tables/term2term.txt",
		root_node
	)	
	
#	stop("Bis hierhin und nicht weiter! Geht die Creation vom Randomset")
	
	# category test			
	FUN2(
		paste(directory, "/randset_out",sep=""),
		paste(directory,"/category_test_out",sep=""),
		1,
		root_node
	)	
}	
	

