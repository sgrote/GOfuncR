
# get genes from genomic regions

# 'genes' is the input-dataframe for go_enrich, in this case it's genomic regions like 9:123-1225
# output: list of
    # classic single genes go_enrich input-dataframe
    # gene_coords chr, start, end, gene for all genes (test+bg)
# side-effect: write regions to file for FUNC
# requires database to be laoded

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
    
    # load database
    if (!silent){
        message(paste("load database '", coord_db,"'...",sep=""))
    }
    if (!suppressPackageStartupMessages(suppressMessages(require(coord_db, character.only=TRUE)))){
        stop(paste0("database '" ,coord_db, "' is not installed. Please install it from bioconductor."))
    }
    if (!silent){
        message(paste("find genes in input-regions using database '", coord_db,"'...",sep=""))
    }
    
    # convert to GRanges
    test_range = GRanges(test_regions[,1], IRanges(test_regions[,2], test_regions[,3]))
    bg_range = GRanges(bg_regions[,1], IRanges(bg_regions[,2], bg_regions[,3]))
    
    # check if OrganismDb or OrgDb/TxDb
    if (coord_db == anno_db){
        gene_identifier = "SYMBOL"
    } else {
        gene_identifier = "GENEID"
    }
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
        stop("Candidate regions do not contain protein-coding genes.")
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
        stop("All values of the 'genes[,2]'-input are 1. Using chromosomal regions as input requires defining background regions with 0.")
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


# taken from https://gist.github.com/mtmorgan/bcacbea1b46445769f1cb91f87e25c30
geneRanges = function(db=Homo.sapiens, column="SYMBOL"){
    g = genes(db, columns=column)
    col = mcols(g)[[column]]
    genes = granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] = as.character(unlist(col))
    return(genes)
}

# find overlaps of genes, ranges and convert to data.frame chr, start, end, gene
get_genes_from_regions = function(gene_coords, ranges){
    genes = subsetByOverlaps(gene_coords, ranges)
    out = data.frame(chr=seqnames(genes), start=start(genes), end=end(genes), gene=elementMetadata(genes)[,1])
    return(out)
}

# convert EntrezID from TxDb to symbol using orgDb
entrez_to_symbol = function(entrez, orgDb=org.Hs.eg.db){
    symbol = suppressMessages(select(orgDb, keys=as.character(entrez), columns=c("ENTREZID","SYMBOL"), keytype="ENTREZID"))
    # just double-check
    if(any(symbol[,1] != entrez)) {stop("Unexpected order in entrez_to_symbol.")}
    return(symbol)
}

# load database
load_db = function(db, silent=FALSE){
    if (!silent){
        message(paste("load database '", db,"'...",sep=""))
    }
    if (!suppressPackageStartupMessages(suppressMessages(require(db, character.only=TRUE)))){
        stop(paste0("database '" ,db, "' is not installed. Please install it from bioconductor."))
    }
}

# take the lowest transcript start and the highest end (cols in: gene, chr, start, end)
combine_tx = function(coords){
    c1 = aggregate(TXSTART ~ SYMBOL + TXCHROM, coords, min)
    c2 = aggregate(TXEND ~ SYMBOL + TXCHROM, coords, max)
    out = merge(c1, c2)[,c(2:4,1)] # move gene to last column
    return(out)
}

# get chr, start, stop, symbol for all genes in coord_db
get_all_coords = function(coord_db="Homo.sapiens", anno_db="Homo.sapiens", silent=FALSE){
    load_db(coord_db, silent)
    if (!silent){
        message(paste("find gene coordinates using database '", coord_db,"'...",sep=""))
    }   
    if (coord_db == anno_db){
        # OrganismDb
        symbols = keys(get(coord_db), keytype="SYMBOL")
        coords = suppressMessages(select(get(coord_db), keys=symbols, columns=c("TXCHROM", "TXSTART", "TXEND", "SYMBOL"), keytype="SYMBOL"))
    } else {
        # OrgDb/TxDb (combine possible different Entrez per Symbol)
        load_db(anno_db, silent)
        entrez = keys(get(coord_db), keytype="GENEID")
        coords = suppressMessages(select(get(coord_db), keys=entrez, columns=c("TXCHROM", "TXSTART", "TXEND", "GENEID"), keytype="GENEID"))
        coords[,1] = entrez_to_symbol(coords[,1], get(anno_db))[,2]
        colnames(coords)[1] = "SYMBOL"
    }
    # remove invalid chroms like "chr6_qbl_hap6"
    coords = coords[grepl("^chr[0-9XYMT]*$", coords[,2]),]
    # maximum transcript range
    coords = combine_tx(coords)
    # remove genes that are still duplicated (on different chroms, mostly X/Y)
    coords = coords[!(coords[,4] %in% coords[duplicated(coords[,4]),4]) ,]
    
    return(coords)
}

# get chr, start, stop, symbol for input symbols
get_gene_coords = function(symbols, coord_db="Homo.sapiens", anno_db="Homo.sapiens", silent=FALSE){
    load_db(coord_db, silent)
    if (coord_db == anno_db){
        # OrganismDb
        coords = suppressMessages(select(get(coord_db), keys=symbols, columns=c("TXCHROM", "TXSTART", "TXEND", "SYMBOL"), keytype="SYMBOL"))
    } else {
        # OrgDb/TxDb (combine possible different Entrez per Symbol)
        load_db(anno_db, silent)
        entrez = suppressMessages(select(get(anno_db), keys=symbols, columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL"))
        en_coords = suppressMessages(select(get(coord_db), keys=entrez[,2], columns=c("TXCHROM", "TXSTART", "TXEND", "GENEID"), keytype="GENEID"))
        coords = merge(entrez, en_coords, by.x="ENTREZID", by.y="GENEID")[,2:5]
    }
    # remove invalid chroms like "chr6_qbl_hap6"
    coords = coords[grepl("^chr[0-9XYMT]*$", coords[,2]),]
    # maximum transcript range
    coords = combine_tx(coords)
    # remove genes that are still duplicated (on different chroms, mostly X/Y)
    coords = coords[!(coords[,4] %in% coords[duplicated(coords[,4]),4]) ,]
    
    return(coords)
}





