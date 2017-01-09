library(FuncBlocks)
set.seed(123)

## run Func
gene_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
go_res = go_enrich(genes, n_randset=10)

## run Func - defined background
candi_gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
bg_gene_ids = c('FGR', 'NPHP1', 'DRD2', 'ABCC10', 'PTBP2', 'JPH4', 'SMARCC2', 'FN1', 'NODAL', 'CYP1A2', 'ACSS1', 'CDHR1', 'SLC25A36', 'LEPR', 'PRPS2', 'TNFAIP3', 'NKX3-1', 'LPAR2', 'PGAM2', 'GAPDHS')
genes = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
names(genes) = c(candi_gene_ids, bg_gene_ids)
go_bg = go_enrich(genes, n_randsets=10)


##### plot_odds ratio
# normal input
plot_odds_ratio(go_res, fwer_threshold=0.02)
plot_odds_ratio(go_bg,fwer_threshold=0.8)
plot_odds_ratio(go_res, go_ids=c('GO:0072025','GO:0072221','GO:0072235', 'GO:0044765'))
plot_odds_ratio(go_bg, go_ids=c('GO:0005634','GO:0004945','0.05309471','GO:0008289','GO:0005737','GO:0071495'))
# erroneous input
plot_odds_ratio(go_bg, go_ids=c('GO:0000009', 'GO:0000010', 'GO:0000014')) # no genes annotated
plot_odds_ratio(go_bg, fwer_threshold=0.001) # no gos with fwer below threshold
plot_odds_ratio("bla")
plot_odds_ratio(go_bg, fwer_threshold=2) # all GOs

##### corner case
# only one GO
plot_odds_ratio(go_bg, go_ids=c('GO:0005623'))
# no candidate annotated
plot_odds_ratio(go_bg, go_ids=c('GO:0000166', 'GO:0000287', 'GO:0000981'))
