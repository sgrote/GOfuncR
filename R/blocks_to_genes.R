
# get genes from genomic regions

# 'genes' is the input-dataframe for go_enrich, in this case it's genomic regions like 9:123-1225
# output: list of
    # classic single genes go_enrich input-dataframe
    # gene_coords chr, start, end, gene for all genes (test+bg)
# side-effect: write regions to file for FUNC

blocks_to_genes = function(directory, genes, anno_db="Homo.sapiens", coord_db="Homo.sapiens", circ_chrom=FALSE, silent=FALSE){

    # check regions are valid, remove unused chroms for circ_chrom,
    # get two bed-dataframes back (candidate and background)
    regions = check_regions(genes, circ_chrom)
    
    # add "chr" substring (like in databases)
    regions = lapply(regions, function(x) {x[,1] = paste0("chr", x[,1]); return(x)})
    test_regions = regions[[1]]
    bg_regions = regions[[2]]
    if (!silent){
        message("Candidate regions:")
        print(test_regions)
        message("Background regions:")
        print(bg_regions)
    }
    
    # load databases and check if OrganismDb or OrgDb/TxDb
    load_db(anno_db, silent)
    if (coord_db == anno_db){
        gene_identifier = "SYMBOL"
    } else {
        gene_identifier = "GENEID"
        load_db(coord_db, silent)
    }
    
    if (!silent){
        message(paste("find genes in input-regions using database '", coord_db,"'...",sep=""))
    }
    
    # convert to GRanges
    test_range = GRanges(test_regions[,1], IRanges(test_regions[,2], test_regions[,3]))
    bg_range = GRanges(bg_regions[,1], IRanges(bg_regions[,2], bg_regions[,3]))
    
    # get overlapping genes
    all_genes = suppressMessages(geneRanges(get(coord_db), column=gene_identifier))
    test_genes = get_genes_from_regions(all_genes, test_range)
    bg_genes = get_genes_from_regions(all_genes, bg_range)
    # convert EntrezID from TxDb to symbol using orgDb
    if (gene_identifier == "GENEID"){
        test_genes[,4] = entrez_to_symbol(test_genes[,4], get(anno_db))[,2]
        bg_genes[,4] = entrez_to_symbol(bg_genes[,4], get(anno_db))[,2]
        test_genes = test_genes[!is.na(test_genes[,4]),]
        bg_genes = bg_genes[!is.na(test_genes[,4]),]
    }
    
    # check that candidate and background contain genes
    if (nrow(test_genes) == 0){
        stop("Candidate regions do not contain any genes.")
    }
    if (nrow(bg_genes) == 0){
        stop("Background regions do not contain any genes.")
    }
    
    ## write candidate and background bed-files for C++
    # merge candidate into background for randomsets
    full_bg_range = suppressWarnings(reduce(c(test_range, bg_range)))
    full_bg_regions = data.frame(chr=seqnames(full_bg_range), start=start(full_bg_range), end=end(full_bg_range))
    # avoid scientific notation in regions (read in c++)
    test_regions = format(test_regions, scientific=FALSE, trim=TRUE)
    full_bg_regions = format(full_bg_regions, scientific=FALSE, trim=TRUE)
    # write regions to files
    write.table(test_regions,file=paste(directory, "_test_regions.bed",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
    write.table(full_bg_regions,file=paste(directory, "_bg_regions.bed",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

    # combine gene coords of background and test regions (to create input for FUNC in go_enrich)
    gene_coords = unique(rbind(test_genes, bg_genes))

    # convert to classic go_enrich input avoiding double-assignment of candidate/background
    genes_df = data.frame(genes=gene_coords$gene, score=rep(0,nrow(gene_coords)))
    genes_df[,1] = as.character(genes_df[,1])
    genes_df[genes_df$genes %in% test_genes[,4], 2] = 1
    
    return(list(genes_df, gene_coords))
}




# input: 
    # dataframe with genomic regions (chr:from-to) and 1/0
    # circ_chrom-option T/F
# output: list with elements
    # test-regions (bed)
    # background_regions (merged with test-regions) (bed)
check_regions = function(genes, circ_chrom){
    
    # check that background region is specified
    if (all(genes[,2]==1)){
        stop("All values of the genes[,2] input are 1. Using chromosomal regions as input requires defining background regions with 0.")
    }
    
    # remove invalid ranges on (mt)-chromosome
    inval = genes[! grepl("^[0-9XYMT]*:[0-9]*-[0-9]*$", genes[,1]),1]
    if (length(inval) > 0){
        stop(paste("Invalid regions:", paste(inval, collapse=", ")))
    }
    
    # convert coordinates from 'genes'-names to bed-format
    genes[,1] = as.character(genes[,1])
    bed = do.call(rbind, strsplit(genes[,1], "[:-]"))
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

    return(list(test_reg, bg_reg))    
}





