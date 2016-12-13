# remove package from R session
detach("package:FuncBlocks", unload = TRUE)
#library.dynam.unload("FuncBlocks", system.file(package = "FuncBlocks"))

library(FuncBlocks)
setwd("/r1/people/steffi_grote/R_packages/FuncBlocks_package")

## get name 
get_GO_names(c("GO:0051082", "GO:123", "GO:0042254", "GO:0000109"))


# TODO: auf neues genes-output anpassen

###### Blocks

background = read.table("Ben_background_regions.bed")
candidate = read.table("Ben_neandertal_deserts.bed")
desmap = rbind(candidate, background)
genes = c(rep(1,nrow(candidate)),rep(0,nrow(background)))
names(genes) = c(paste(desmap[,1],":",desmap[,2],"-",desmap[,3],sep=""))


# normal run
res1 = go_enrich(genes)
res2 = go_enrich(genes, circ_chrom=TRUE)
head(res1$results)
head(res2$results)
head(res1$genes)
head(res2$genes)

# check that genes from regions yield same p-val when passed directly to go_enrich
res11 = go_enrich(res1$genes)
res22 = go_enrich(res2$genes)
head(res1$results)
head(res11$results)
head(res2$results)
head(res22$results)

## Checks, Warnings and Error messages

# wilcoxon -> does not work with regions - could that make sense? probably not -> too many ties?
willi = as.integer(runif(length(genes),1,300))
names(willi) = names(genes)
res = go_enrich(willi, test="wilcoxon", n_randsets=5)

# background region too small that non-overlapping random block placement might fail
tight = genes[c(1:4,36)]
res = go_enrich(tight,n_randsets=5)

# background region von vornherein too small
too_small = genes[1:5]
res = go_enrich(too_small)
# test circ_chrom error, no background for candidate on chrom
res = go_enrich(too_small, circ_chrom=TRUE)

# test circ_chrom error, background < candidate on chrom
too_small = c(rep(1,4),rep(0,5)) 
names(too_small) = c("X:0-2","2:0-3","2:5-10","13:0-20",  "2:5-10","13:0-10","X:0-1","4:0-10","2:10-12") 
res = go_enrich(too_small,n_randsets=5, circ_chrom=TRUE)

# candidate regions overlap
over = genes
names(over)[3] = "8:60500000-62400000" # one region inside another
res = go_enrich(over,n_randsets=5)

# background regions overlap
over = genes
names(over)[7] = "1:27000000-161100000" # normal overlap
names(over)[32] = "3:49700000-94600000" # normal overlap
res = go_enrich(over)
res = go_enrich(over,circ_chrom=TRUE)

# start > stop
err_genes = genes
names(err_genes)[3] = "7:124700000-113600000"
names(err_genes)[10] = "1:13200000-1500000"
res = go_enrich(err_genes)

# test error when no background region specified
err_genes = err_genes[err_genes == 1]
res = go_enrich(err_genes)

## chr21 as chr in region input - leads to "No requested genes in data. "
err_genes = genes
names(err_genes)[1] = "chr1:104000000-114900000" 
res = go_enrich(err_genes)

## warning if gene_len is TRUE but not used
res = go_enrich(genes, n_randsets=5, gene_len=TRUE)




###### non-Blocks
detach("package:FuncBlocks", unload = TRUE)
library(FuncBlocks)
setwd("/r1/people/steffi_grote/R_packages/FuncBlocks_package")

test_genes = paste(rep('FABP', 5), c(2,4:7), sep='')
test_genes = c(test_genes, 'HSPB2' ,'LINC00239', 'TESTI') # die sollten keine coordinates haben
bg_genes = c('NCAPG', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'MIRLET7BHG')
genes = c(rep(1,length(test_genes)), rep(0,length(bg_genes)))
names(genes) = c(test_genes, bg_genes)

res3=go_enrich(genes, n_randsets=5)
res4=go_enrich(genes, n_randsets=5, gene_len=TRUE)
print(head(res3[[1]]))
print(head(res4[[1]]))

# no background
nb_genes = genes[genes==1]
res5=go_enrich(nb_genes, n_randsets=5)

# test warning if circ_chrom is TRUE but not used
res = go_enrich(genes, n_randsets=5, circ_chrom=TRUE)

# wilcoxon
willi = as.integer(runif(length(genes),1,300))
names(willi) = names(genes)
res = go_enrich(willi, test="wilcoxon", n_randsets=5)
#res = go_enrich(willi, test="wilcoxon", gene_len=TRUE, n_randsets=5)

# wilcoxon ties
willi = c(1,3,4,3,8,10,10,10,20,20) # last gene not annotated to last root - no real tie
names(willi) = c("AGTR1","ANO1","BTBD3","C21orf59","CACNG2","FABP2","FABP4","FABP5","FABP6","HSPB2")
res = go_enrich(willi, test="wilcoxon", n_randsets=5)


randout = read.table("./tmp/randset_out", skip=1, header=T)
colnames(randout) = gsub("\\.",":", colnames(randout))
rowSums(randout)
table(unlist(c(randout))) # max(sum-scores) = sum of all scores 1-9 (9 anno genes -> 9*5 = 45)
root_node_ids= c("GO:0003674","GO:0008150","GO:0005575")
root = root_node_ids[root_node_ids %in% colnames(randout)]
randout[,root]












