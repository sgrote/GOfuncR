
## take a bed-like-file and merge overlapping regions
## (for merging candidate into background regions in ABAEnrichment)


# merging for one chrom, input is sorted
merge_chrom = function(chrbed){	
	# go through regions and merge with neighbor region
	i = 1
	while (i < nrow(chrbed)){ 	# last line cannot be changed
		# (a) no overlap    (stop before next start)
		if (chrbed[i,3] < chrbed[i+1,2]){
			i = i+1
		} 
		# (b) including the next region
		else if (chrbed[i,3] >= chrbed[i+1,3]){
			chrbed = chrbed[-(i+1),]
		}		
		# (c) included in next region (sorted by start,end -> must be included if starts are identical)
		else if (chrbed[i,2] == chrbed[i+1,2]){
			chrbed = chrbed[-i,]
		}		
		# (d) overlapping next region
		else if (chrbed[i,3] >= chrbed[i+1,2]){
			chrbed[i,3] = chrbed[i+1,3]
			chrbed = chrbed[-(i+1),]
		}
	}
	return(chrbed)
}

merge_bed = function(bed){
	# check that start < stop (athough this is already done in get_genes_from_regions...)
	reverse_indi = bed[,2] > bed[,3]
	if (sum(reverse_indi) > 0){
		reverse = paste(names(genes)[reverse_indi], collapse=", ")
		stop(paste("Invalid regions: ", reverse, ".\n  In 'chr:start-stop' start < stop is required.", sep=""))
	}
	# sort
	bed = bed[order(bed[,1],bed[,2],bed[,3]),]
	# merge
	out = do.call(rbind,by(bed, bed[,1], merge_chrom))
	rownames(out) = 1:nrow(out)
	return(out)
}


#### for TESTING
### print before and after merging
#plot_bed = function(bed){
#	plot(1,type="n",ylim=range(chroms),xlim=range(c(starts,stops)),xlab="",ylab="")
#	apply(bed,1,function(reg) lines(c(reg[2],reg[3]),c(reg[1],reg[1])+runif(1,-0.1,0.1)))
#}
#chroms = c(1,1,1,2,2,2,2,2,2,3,3,4,4,4,5,5,5,5)
#starts = c(2,4,10,2,2,2,7,7,7,1,3,8,2,5,1,7,2,3)
#stops = c(3,10,12,5,5,5,9,8,10,4,5,10,7,6,2,9,4,5)
#bed = data.frame(chroms,starts,stops)
#par(mfrow=c(1,2))
#plot_bed(bed)
#bed_merged = merge_bed(bed)
#plot_bed(bed_merged)



