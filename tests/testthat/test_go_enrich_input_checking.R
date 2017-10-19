
## context("main function go_enrich input checking")

# valid genes input
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
    'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2', 'QUATSCH')
is_candidate = rep(1, length(gene_ids))
input_hyper = data.frame(gene_ids, is_candidate, stringsAsFactors=FALSE)
# multi-assign
input_multi = rbind(input_hyper, input_hyper[1:2,])
input_multi[1:2,2] = 0
# not 0/1
input_non_bin = input_hyper 
input_non_bin[3,2] = 2

## general input
test_that("general genes input gets checked",{
	expect_error(go_enrich(1:10),"Please provide a data frame as 'genes' input")
	expect_error(go_enrich(input_hyper, n_randset="abc"),
	    "Please define 'n_randsets' as a positive integer.")
	expect_error(go_enrich(input_hyper, circ_chrom=TRUE),
		"Argument 'circ_chrom = TRUE' can only be used with 'regions = TRUE'.")
	expect_error(go_enrich(input_hyper, test="invalid_test"),
		"Not a valid test. Please use 'hyper', 'wilcoxon', 'binomial' or 'contingency'.")
	expect_error(go_enrich(input_multi),
		"Genes with multiple assignment in input: APOL4, NCAPG")
	expect_error(go_enrich(input_non_bin),
		"Please provide only 1/0-values in 2nd column of 'genes'-input for hypergeometric test.")
})

## hyper
# one gene without annotation
genes_no_anno = data.frame(gene_ids='QUATSCH1' , scores=1) 
candi_no_anno = data.frame(gene_ids=c('QUATSCH1','APOL4') , scores=c(1,0))
bg_no_anno = data.frame(gene_ids=c('QUATSCH1','APOL4') , scores=c(0,1))

test_that("hyper genes input gets checked",{
	expect_error(go_enrich(genes_no_anno), "No GO-annotations for input genes.")
	expect_error(go_enrich(candi_no_anno), "No requested candidate genes in data.")
	expect_error(go_enrich(bg_no_anno), "No requested background genes in data.")
})

# regions
# background < candidate on chrom, in blocks, overlapping input regions
too_small = c(rep(1,4),rep(0,5)) 
names(too_small) = c("X:0-2","2:0-3","2:5-10","13:0-20",  "2:5-15","13:0-10","X:0-3","4:0-10","2:10-12") 
overlap = too_small
names(overlap) = c("X:0-1","2:0-3","2:5-10","13:0-20",  "2:5-20","13:0-30","X:0-3","4:0-10","2:10-12")
overlap2 = too_small
names(overlap2)[2] = "2:0-8"
tight = c(rep(1,4),0)
names(tight) = c("1:104000000-114900000", "3:76500000-90500000", "7:113600000-124700000", "8:54500000-65400000", "5:0-4700000")
no_bg = too_small[1:4]
reverse = too_small
names(reverse)[3] = "2:15-10"
no_can_genes = c(1, rep(0,6))
names(no_can_genes) = c("1:10-20", "7:1300000-56800000", "7:74900000-148700000", "8:7400000-44300000", "8:47600000-146300000", "9:0-39200000", "9:69700000-140200000")
no_bg_genes = c(1,0)
names(no_bg_genes) = c("8:82000000-83000000", "21:1-3000000")

test_that("input_regions are checked - blocks",{
    expect_that(go_enrich(overlap, regions=TRUE), throws_error("Background regions overlap: 2:5-20, 2:10-12"))
    expect_error(go_enrich(overlap2, regions=TRUE), "Candidate regions overlap: 2:0-8, 2:5-10")
    expect_error(go_enrich(tight, regions=TRUE), "Background regions too small.")
    expect_error(go_enrich(no_bg, regions=TRUE), "All values of the genes") # longer message fails ("["?)
	expect_error(go_enrich(reverse, regions=TRUE), 
	    "Invalid regions: 2:15-10.\n  In 'chr:start-stop' start < stop is required.")
	expect_error(go_enrich(no_can_genes, regions=TRUE), "Candidate regions do not contain any genes.")
	expect_error(go_enrich(no_bg_genes, regions=TRUE), "Background regions do not contain any genes.")
})    

test_that("input_regions are checked - circ_chrom",{
    expect_error(go_enrich(overlap, circ_chrom=TRUE, regions=TRUE),
        "Background regions overlap: 2:5-20, 2:10-12")
    expect_error(go_enrich(overlap2, circ_chrom=TRUE, regions=TRUE),
        "Candidate regions overlap: 2:0-8, 2:5-10")
    expect_error(go_enrich(no_bg, circ_chrom=TRUE, regions=TRUE), 
        "All values of the genes")
    expect_error(go_enrich(reverse, circ_chrom=TRUE, regions=TRUE), 
        "Invalid regions: 2:15-10.\n  In 'chr:start-stop' start < stop is required.")
})

# wilcox
high_score_genes = c('GCK', 'QUATSCH')
low_score_genes = c('CACNG2')
gene_scores = 1:3
input_willi = data.frame(gene_ids = c(high_score_genes, low_score_genes), gene_scores)
input_extra_col = cbind(input_willi, input_willi[,2])

test_that("willi genes input gets checked",{
	expect_error(go_enrich(input_willi[1:2,], test="wilcoxon"),
	    "Less than 2 genes have annotated GO-categories.")
	expect_error(go_enrich(input_extra_col, test="wilcoxon"),
	    "Please provide a data frame with 2 columns.") # again [ character can not be checked

})

# binomial
gene_ids = c('G6PD', 'GCK', 'CACNG2', 'QUATSCH')
A_counts = 1:4
B_counts = 2:5
input_binom = data.frame(gene_ids, A_counts, B_counts)

test_that("binom genes input gets checked",{
	expect_error(go_enrich(input_binom[4,], test="binomial"),
	    "None of the genes entered are present in the SYMBOL column of 'Homo.sapiens'. Check head")
	expect_error(go_enrich(input_binom[,c(1:3,3)], test="binomial"),
	    "Please provide a data frame with columns ") # again [ character can not be checked
})

# contingency
gene_ids = c('G6PD', 'GCK', 'GYS1', 'QUATSCH')
A_counts = 1:4
B_counts = 2:5
C_counts = 1:4
D_counts = 2:5
input_conti = data.frame(gene_ids, A_counts, B_counts, C_counts, D_counts)

test_that("conti genes input gets checked",{
	expect_error(go_enrich(input_conti[4,], test="contingency"),
	    "None of the genes entered are present in the SYMBOL column of 'Homo.sapiens'. Check head")
	expect_error(go_enrich(input_conti[,c(1:3,3)], test="contingency"),
	    "Please provide a data frame with columns ") # again [ character can not be checked
})

