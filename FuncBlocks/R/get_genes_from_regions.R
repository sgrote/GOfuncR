
# input: 
	# 0/1-vector with genomic regions as names (chr:from-to)
	# gene_coords (bed)
	# circ_chrom-option T/F
# output: list with elements
	# test-regions (bed)
	# background_regions (bed)
	# genes-vector ("normal hyper-input" for genes from test-regions) 

get_genes_from_regions = function(genes, gene_pos, circ_chrom){
	
	# convert coordinates from 'genes'-names to bed-format 
	bed = do.call(rbind, strsplit(names(genes), "[:-]"))
#	bed[,1] = substring(bed[,1], 4)  ## if chrom is given as chr21 instead of 21 
	bed = as.data.frame(bed)
	bed[,2:3] = apply(bed[,2:3], 2, as.numeric)
	bed[,1] = as.character(bed[,1])
	
	# check that start < stop
	reverse_indi = bed[,2] > bed[,3]
	if(sum(reverse_indi) > 0){
		reverse = paste(names(genes)[reverse_indi], collapse=", ")
		stop(paste("Invalid regions: ", reverse, ".\n  In 'chr:start-stop' start < stop is required.", sep=""))
	}
		
	# split in test and background
	test_reg = bed[genes==1,]
	bg_reg = bed[genes==0,]	
	
	## test that regions are non-overlapping (separately for candidate and background)
	# candidate
	overlap_indis = c()
	for (i in 1:nrow(test_reg)){ 
		# if (chrom=chrom & (start inside | end inside | including)) any other region
		if (any(test_reg[i,1] == test_reg[,1] & ((test_reg[i,2] > test_reg[,2] & test_reg[i,2] < test_reg[,3]) | (test_reg[i,3] > test_reg[,2] & test_reg[i,3] < test_reg[,3]) | (test_reg[i,2] < test_reg[,2] & test_reg[i,3] > test_reg[,3])))){
			overlap_indis = c(overlap_indis, i)
		}	
	}
	if (length(overlap_indis) > 0){
		over = paste(names(genes[genes==1])[overlap_indis], collapse=", ")
		stop(paste("Candidate regions overlap: ", over, sep=""))
	}
	# background
	overlap_indis = c()
	for (i in 1:nrow(bg_reg)){ 
		if (any(bg_reg[i,1] == bg_reg[,1] & ((bg_reg[i,2] > bg_reg[,2] & bg_reg[i,2] < bg_reg[,3]) | (bg_reg[i,3] > bg_reg[,2] & bg_reg[i,3] < bg_reg[,3]) | (bg_reg[i,2] < bg_reg[,2] & bg_reg[i,3] > bg_reg[,3])))){
			overlap_indis = c(overlap_indis, i)
		}	
	}
	if (length(overlap_indis) > 0){
		over = paste(names(genes[genes==0])[overlap_indis], collapse=", ")
		stop(paste("Background regions overlap: ", over, sep=""))
	}
		
	# sort  (mixedorder does not work with multiple columns) 
#	test_reg = test_reg[order(test_reg[,1], test_reg[,2]),] 
#	bg_reg = bg_reg[order(bg_reg[,1], bg_reg[,2]),] 
	test_reg = test_reg[order(test_reg[,2]),] 
	bg_reg = bg_reg[order(bg_reg[,2]),] 
	test_reg = test_reg[mixedorder(test_reg[,1]),] 
	bg_reg = bg_reg[mixedorder(bg_reg[,1]),] 
	
	# if rolling chrom: remove unused bg chroms and warn, check that all candidate chroms have bg
	if(circ_chrom == TRUE){				
		if(!(all(bg_reg[,1] %in% test_reg[,1]))){
			not_used = unique(bg_reg[!(bg_reg[,1] %in% test_reg[,1]),1])
			not_used = paste(not_used, collapse=", ")
			warning(paste("Unused chromosomes in background regions: ", not_used, ".\n  With circ_chrom=TRUE only background regions on the same chromosome as a candidate region are used.",sep=""))
			bg_reg = bg_reg[bg_reg[,1] %in% test_reg[,1],]
		}			
		if(!(all(test_reg[,1] %in% bg_reg[,1]))){
			wo_bg = unique(test_reg[!(test_reg[,1] %in% bg_reg[,1]),1])
			wo_bg = paste(wo_bg, collapse=", ")
			stop(paste("No background region for chromosomes: ",  wo_bg, ".\n  With circ_chrom=TRUE only background regions on the same chromosome as a candidate region are used." ,sep=""))
		}
		# check that background is big enough on all chroms, bg and test chroms are the same and have the same order 
		bg_length = tapply(bg_reg[,3] - bg_reg[,2], bg_reg[,1], sum)
		test_length = tapply(test_reg[,3] - test_reg[,2], test_reg[,1], sum)
		too_big_indi = test_length > bg_length
		if (sum(too_big_indi) > 0){
			too_big = names(bg_length)[too_big_indi]
			too_big = paste(mixedsort(too_big), collapse=", ")
			stop(paste( "Sum of candidate regions is bigger than sum of background regions on chromosomes: ", too_big, sep=""))
		}
				
	} else {  # normal blocks option
		# check that biggest test_region is not bigger than biggest bg_region
		max_bg = max(bg_reg[,3] - bg_reg[,2])
		too_big_indi = (test_reg[,3] - test_reg[,2]) > max_bg
		if (sum(too_big_indi) > 0){
			too_big = paste(paste(test_reg[,1],":",test_reg[,2],"-",test_reg[,3],sep="")[too_big_indi], collapse=", ")
			stop(paste( "Candidate regions bigger than any background region:\n  ", too_big, sep=""))
		}
		# sort candidate regions by length (better chances that random placement works with small bg-regions)
		test_reg = test_reg[order(test_reg[,3] - test_reg[,2], decreasing=T),]
	}
	
	# get genes overlapping background-regions
	bg_genes = c()
	for (i in 1:nrow(bg_reg)){
		bg_genes = c(bg_genes, gene_pos[gene_pos[,"chr"]==bg_reg[i,1] & ((gene_pos[,"start"] >= bg_reg[i,2] & gene_pos[,"start"] < bg_reg[i,3]) | (gene_pos[,"end"] >= bg_reg[i,2] & gene_pos[,"end"] < bg_reg[i,3]) |  (gene_pos[,"start"] <= bg_reg[i,2] & gene_pos[,"end"] >= bg_reg[i,3])), "hgnc_symbol"])		
	}
	# check that bg-region contains genes
	# (if no bg-genes here, all non-candidate genes would be background in go_enrich -> unwanted)
	if (length(bg_genes)==0){
		stop("Background regions do not contain protein-coding genes.")
	}
	
	# get genes overlapping test-regions
	test_genes = c()
	for (i in 1:nrow(test_reg)){
		test_genes = c(test_genes, gene_pos[gene_pos[,"chr"]==test_reg[i,1] & ((gene_pos[,"start"] >= test_reg[i,2] & gene_pos[,"start"] < test_reg[i,3]) | (gene_pos[,"end"] >= test_reg[i,2] & gene_pos[,"end"] < test_reg[i,3]) | (gene_pos[,"start"] <= test_reg[i,2] & gene_pos[,"end"] >= test_reg[i,3])), "hgnc_symbol"])		
	}
	# check that test-region contains genes
	if (length(test_genes)==0){
		stop("Candidate regions do not contain protein-coding genes.")
	}

	# convert to classic "genes" func-input-vector 
	gene_names = unique(c(test_genes, bg_genes))
	genes_vec = rep(0, length(gene_names))
	names(genes_vec) = gene_names
	genes_vec[names(genes_vec) %in% test_genes] = 1
	
	return(list(test_reg, bg_reg, genes_vec))	
}
	
	
	
	
	
	
	
	
