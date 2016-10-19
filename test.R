
library(FuncBlocks)
setwd("/r1/people/steffi_grote/R_packages/FuncBlocks_package")

#### Blocks
background = read.table("Ben_background_regions.bed")
candidate = read.table("Ben_neandertal_candidate.bed")
desmap = rbind(candidate, background)
genes = c(rep(1,nrow(candidate)),rep(0,nrow(background)))
names(genes) = c(paste(desmap[,1],":",desmap[,2],"-",desmap[,3],sep=""))

res1 = go_enrich(genes,cutoff_quantiles=c(0.2,0.7),n_randsets=5)
res2 = go_enrich(genes,cutoff_quantiles=c(0.2,0.7),n_randsets=5, circ_chrom=TRUE)
res1 = go_enrich(genes,dataset='5_stages',cutoff_quantiles=c(0.2,0.7),n_randsets=5)
res2 = go_enrich(genes,dataset='5_stages',cutoff_quantiles=c(0.2,0.7),n_randsets=5, circ_chrom=TRUE)

# normal run
res1 = go_enrich(genes)
res2 = go_enrich(genes, circ_chrom=TRUE)
print(head(res1[[1]]))
print(head(res2[[1]]))

# background region too small that non-overlapping random block placement might fail
tight = genes[c(1:4,36)]
res1 = go_enrich(tight,cutoff_quantiles=c(0.2,0.7),n_randsets=5)

# background region von vornherein too small
too_small = genes[1:5]
res1 = go_enrich(too_small)
# test circ_chrom error, no background for candidate on chrom
res1 = go_enrich(too_small, circ_chrom=TRUE)

# test circ_chrom error, background < candidate on chrom
too_small = c(rep(1,4),rep(0,5)) 
names(too_small) = c("X:0-2","2:0-3","2:5-10","13:0-20",  "2:5-10","13:0-10","X:0-1","4:0-10","2:10-12") 
res1 = go_enrich(too_small,cutoff_quantiles=c(0.2,0.7),n_randsets=5, circ_chrom=TRUE)

# candidate regions overlap
over = genes
names(over)[3] = "8:60500000-62400000" # one region inside another
res1 = go_enrich(over,cutoff_quantiles=c(0.2,0.7),n_randsets=5)

# background regions overlap
over = genes
names(over)[7] = "1:27000000-161100000" # normal overlap
names(over)[32] = "3:49700000-94600000" # normal overlap
res2 = go_enrich(over)
res2 = go_enrich(over,circ_chrom=TRUE)

# start > stop
err_genes = genes
names(err_genes)[3] = "7:124700000-113600000"
names(err_genes)[10] = "1:13200000-1500000"
res1 = go_enrich(err_genes)

# test error when no background region specified
err_genes = err_genes[err_genes == 1]
res1 = go_enrich(err_genes)

## chr21 as chr in region input - leads to "No requested genes in data. "
err_genes = genes
names(err_genes)[1] = "chr1:104000000-114900000" 
res1 = go_enrich(err_genes)

## warning if gene_len is TRUE but not used
res1 = go_enrich(genes, cutoff_quantiles=c(0.2,0.7), n_randsets=5, gene_len=TRUE)




#### non-Blocks
test_genes = paste(rep('FABP', 5), c(2,4:7), sep='')
test_genes = c(test_genes, 'HSPB2' ,'LINC00239', 'TESTI') # die sollten kein coordinates haben
bg_genes = c('NCAPG', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'MIRLET7BHG')
genes = c(rep(1,length(test_genes)), rep(0,length(bg_genes)))
names(genes) = c(test_genes, bg_genes)

# normal run, also test warning if gene_len or circ_chrom are TRUE but not used
res = go_enrich(genes, cutoff_quantiles=c(0.2,0.7), n_randsets=5, gene_len=TRUE, circ_chrom=TRUE)


res3=go_enrich(genes, gene_len=TRUE,cutoff_quantiles=c(0.2,0.7),n_randsets=1 )
res4=go_enrich(genes, cutoff_quantiles=c(0.2,0.7),n_randsets=5)
print(head(res3[[1]]))
print(head(res4[[1]]))

# check that quantiles are re-sorted
res3=go_enrich(genes,cutoff_quantiles=c(0.9999, 0.2),n_randsets=5)
res4=go_enrich(genes,cutoff_quantiles=c(0.9999),n_randsets=5)

