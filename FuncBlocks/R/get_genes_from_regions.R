
# input: 
	# 0/1-vector with genomic regions as names (chr:from-to)
	# gene_coords (bed)
	# circ_chrom-option T/F
# output: list with elements
	# test-regions (bed)
	# background_regions (merged with test-regions) (bed)
	# genes-vector ("normal hyper-input" for genes from test-regions) 

get_genes_from_regions = function(genes, gene_pos, circ_chrom){
	
	# convert coordinates from 'genes'-names to bed-format 
	bed = do.call(rbind, strsplit(genes[,1], "[:-]"))
#	bed[,1] = substring(bed[,1], 4)  ## if chrom is given as chr21 instead of 21 
	bed = as.data.frame(bed)
	bed[,2:3] = apply(bed[,2:3], 2, as.numeric)
	bed[,1] = as.character(bed[,1])
	
	# check that start < stop
	reverse_indi = bed[,2] > bed[,3]
	if(sum(reverse_indi) > 0){
		reverse = paste(genes[,1][reverse_indi], collapse=", ")
		stop(paste("Invalid regions: ", reverse, ".\n  In 'chr:start-stop' start < stop is required.", sep=""))
	}
		
	# split in test and background
	test_reg = bed[genes[,2]==1,]
	bg_reg = bed[genes[,2]==0,]	
	
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
		over = paste(genes[genes[,2]==1,1][overlap_indis], collapse=", ")
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
		over = paste(genes[genes[,2]==0,1][overlap_indis], collapse=", ")
		stop(paste("Background regions overlap: ", over, sep=""))
	}
		
	# sort  (mixedorder does not work with multiple columns) 
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
	} else {  # normal blocks option
		# sort candidate regions by length (better chances that random placement works with small bg-regions)
		test_reg = test_reg[order(test_reg[,3] - test_reg[,2], decreasing=TRUE),]
	}
	
	# get genes overlapping background-regions
	bg_genes = c()
	for (i in 1:nrow(bg_reg)){
		bg_genes = c(bg_genes, gene_pos[gene_pos[,1]==bg_reg[i,1] & ((gene_pos[,2] >= bg_reg[i,2] & gene_pos[,2] < bg_reg[i,3]) | (gene_pos[,3] >= bg_reg[i,2] & gene_pos[,3] < bg_reg[i,3]) |  (gene_pos[,2] <= bg_reg[i,2] & gene_pos[,3] >= bg_reg[i,3])), 4])		
	}
	# check that bg-region contains genes
	# (if no bg-genes here, all non-candidate genes would be background in go_enrich -> unwanted)
	if (length(bg_genes)==0){
		stop("Background regions do not contain protein-coding genes.")
	}
	
	# get genes overlapping test-regions
	test_genes = c()
	for (i in 1:nrow(test_reg)){
		test_genes = c(test_genes, gene_pos[gene_pos[,1]==test_reg[i,1] & ((gene_pos[,2] >= test_reg[i,2] & gene_pos[,2] < test_reg[i,3]) | (gene_pos[,3] >= test_reg[i,2] & gene_pos[,3] < test_reg[i,3]) | (gene_pos[,2] <= test_reg[i,2] & gene_pos[,3] >= test_reg[i,3])), 4])
	}
	# check that test-region contains genes
	if (length(test_genes)==0){
		stop("Candidate regions do not contain protein-coding genes.")
	}

	# convert to classic "genes" func-input-dataframe
	gene_names = c(test_genes, bg_genes)
	scores = c(rep(1,length(test_genes)), rep(0,length(bg_genes)))
	genes_df = unique(data.frame(gene_names, scores))

	# NEW: merge candidate into background regions (single-genes-FUNC also implicitly integrates candidate into background genes to choose from in randomsets)
	bg_reg = merge_bed(rbind(bg_reg,test_reg))
	
	return(list(test_reg, bg_reg, genes_df))	
}
	
	
	
	
	
	
	
	
