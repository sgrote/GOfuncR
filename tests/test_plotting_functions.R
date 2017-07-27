
library(FuncBlocks)


################## (1) plot_odds_ratio

set.seed(123)

## run Func
gene_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 
 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
go_res = go_enrich(genes, n_randset=100)

## run Func - defined background
candi_gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
bg_gene_ids = c('FGR', 'NPHP1', 'DRD2', 'ABCC10', 'PTBP2', 'JPH4', 'SMARCC2', 'FN1', 'NODAL', 'CYP1A2', 'ACSS1', 'CDHR1', 'SLC25A36', 'LEPR', 'PRPS2', 'TNFAIP3', 'NKX3-1', 'LPAR2', 'PGAM2', 'GAPDHS')
genes = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
names(genes) = c(candi_gene_ids, bg_gene_ids)
go_bg = go_enrich(genes, n_randsets=100)
## skip one root
set.seed(123)
go_bg_skip = go_enrich(genes, n_randsets=100, domains=c('molecular_function', 'cellular_component'))
go_bg_skip2 = go_enrich(genes, n_randsets=100, domains=c('molecular_function', 'biological_process'))
go_bg_skip3 = go_enrich(genes, n_randsets=100, domains=c('cellular_component'))

# normal input
x = plot_odds_ratio(go_res, fwer_threshold=0.02)
x
plot_odds_ratio(go_bg,fwer_threshold=0.8)
plot_odds_ratio(go_res, go_ids=c('GO:0072025','GO:0072221','GO:0072235', 'GO:0044765')) # check order-preserve
plot_odds_ratio(go_bg, go_ids=c('GO:0005634','GO:0004945','GO:0008289','GO:0005737','GO:0071495'))
# erroneous input
plot_odds_ratio(go_bg, go_ids=c('GO:0000009', 'GO:0000010', 'GO:0000014')) # no genes annotated
plot_odds_ratio(go_bg, fwer_threshold=0.001) # no gos with fwer below threshold
plot_odds_ratio("bla")
plot_odds_ratio(go_res, fwer_threshold=-1)
plot_odds_ratio(go_res, fwer_threshold='bla')
# corner case
plot_odds_ratio(go_bg, fwer_threshold=2) # all GOs
plot_odds_ratio(go_bg, go_ids=c('GO:0005623')) # only one GO
plot_odds_ratio(go_bg, go_ids=c('GO:0000166', 'GO:0000287', 'GO:0000981')) # no candidate annotated
# skipped root node
plot_odds_ratio(go_bg_skip, fwer_threshold=0.9) # no biol_process
plot_odds_ratio(go_bg_skip2, fwer_threshold=0.9) # no cell_comp
plot_odds_ratio(go_bg_skip3, fwer_threshold=0.9) # no cell_comp


################## (2) plot_scores
library(FuncBlocks)
set.seed(123)
high_score_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_score_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = c(sample(20:30, length(high_score_genes)), sample(5:15, length(low_score_genes)))
names(genes) = c(high_score_genes, low_score_genes)
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)
head(go_willi[[1]])
go_willi[[2]]
go_willi_skip = go_enrich(genes, test='wilcoxon', n_randsets=100, domains=c('cellular_component'))

# normal input
w = plot_scores(go_willi, fwer_threshold=0.21)
w
plot_scores(go_willi, go_ids=c('GO:0072025','GO:0072221','GO:0072235', 'GO:0044765')) # check order-preserve
plot_scores(go_willi, go_ids=c('GO:0005634','GO:0004945','GO:0008289','GO:0005737','GO:0071495'))
# corner cases
plot_scores(go_willi, go_ids=c('GO:0005634')) # only one GO
plot_scores(go_willi, fwer_threshold=2) # all GOs
# erroneous input
plot_scores(go_willi, go_ids=c('GO:0000009', 'GO:0000010', 'GO:0000014')) # no genes annotated
plot_scores(go_willi, fwer_threshold=0.001) # no gos with fwer below threshold
plot_scores("bla")
plot_scores(go_willi, fwer_threshold=-1)
plot_scores(go_willi, fwer_threshold='bla')
# skipped root node
plot_scores(go_willi_skip, fwer_threshold=0.9)

# mouse
set.seed(123)
high_score_genes = c('G6pd', 'Gck', 'Gys1', 'Hk2', 'Pygl', 'Slc2a8', 'Ugp2', 'Zwint', 'Engase')
low_score_genes = c('Cacng2', 'Agtr1', 'Ano1', 'Btbd3', 'Mtus1', 'Calb1', 'Gyg1', 'Pax2')
genes = c(sample(20:30, length(high_score_genes)), sample(5:15, length(low_score_genes)))
names(genes) = c(high_score_genes, low_score_genes)
go_willi_mus = go_enrich(genes, test='wilcoxon', n_randsets=100, ref_genome='grcm38')
plot_scores(go_willi_mus, fwer_threshold=0.21)






