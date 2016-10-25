# detect gene identifier given a gene_id of type either entrez, ensembl or hgnc-symbol	
# new: or genomic regions given as chr:from-to

detect_identifier=function(gene_id){
	# de-factor
	if (is.factor(gene_id)){
		gene_id=as.character(gene_id)
	}
	# a) entrez-IDs are numeric
#		numeri=!is.na(suppressWarnings(as.numeric(gene_id)))
#	if (numeri){
#		return("entrezgene")
#	}	
#	# b) ensembl-IDs start with "ENSG" and have 15 characters  	
#	if (nchar(gene_id)==15 && substring(gene_id,1,4)=="ENSG"){
#		return("ensembl_gene_id")
#	}
	# c) genomic region is chr:from-to (TODO: maybe allow also other "chromosomes"?)
	if (grepl("^[0-9XY]*:[0-9]*-[0-9]*$", gene_id)){
		return("blocks")
	}
	# d) neither of the above cases in present in HGNC-symbols
	return("hgnc_symbol")			
}	
