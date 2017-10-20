
# load database
load_db = function(db, silent=FALSE){
    if (!silent){
        message(paste("load database '", db,"'...",sep=""))
    }
    if (!suppressPackageStartupMessages(suppressMessages(require(db, character.only=TRUE)))){
        stop(paste0("database '" ,db, "' is not installed. Please install it from bioconductor."))
    }
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

