\name{plot_odds_ratio}
\alias{plot_odds_ratio}
\title{
Plot proportions of candidate and background genes in enriched or user-defined GO-categories
}
\description{
Uses an object returned from \code{\link{go_enrich}} and extracts the GO-categories that are significantly (given a FWER-threshold) enriched in candidate genes. Alternatively GO-IDs can be defined explicitly.\cr
\code{plot_odds_ratio} then plots for each GO-category the proportion of annotated candidate and background genes in pie charts, whose size is dependent on the total number annotated genes.\cr
The plot also shows pie charts with the proportions of candidate and background genes annotated to the three root nodes (domains) of the gene ontology; as well as the odds ratios for the proportion of candidate genes in the GO-categories given the proportion in the root node.

}

\usage{
plot_odds_ratio(res, fwer_threshold = 0.05, go_ids = NULL)
}
\arguments{
  \item{res}{
an object returned from \code{\link{go_enrich}} (list of 3 elements).
}
  \item{fwer_threshold}{
numeric defining the FWER-threshold. \code{plot_odds_ratio} extracts all GO-categories from \code{res} that have a FWER < \code{fwer_threshold} and plots the proportions of annotated candidate and background genes for those GO-categories. Only used when \code{go_ids} are not specified.
}
  \item{go_ids}{
vector of GO-IDs, e.g. c('GO:0051082', 'GO:0042254'). If defined, this specifies the GO-categories directly and the \code{fwer_threshold} is not used.
}
}

\details{
The odds ratios and the 95%-CIs are computed with Fisher's exact test (two-sided).
}

\value{
A data frame with the results of Fisher's excact test for every GO-category (invisible).
}

\author{
Steffi Grote
}


\seealso{
\code{\link{fisher.test}}\cr
\code{\link{go_enrich}}\cr
\code{\link{get_anno_genes}}\cr
\code{\link{get_names}}\cr
}

\examples{ 
#### Note that argument 'n_randsets' is reduced 
#### to lower computational time in the examples.

## perform an enrichment analysis for some candidate genes
## and plot proportion of candidate and background genes 
## for GO-categories that have a FWER < 0.02
set.seed(123)
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2',
 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
go_res = go_enrich(genes, n_randset=100, silent=TRUE)
plot_odds_ratio(go_res, fwer_threshold=0.02)


## perform an enrichment analysis for some candidate genes
## with defined background
## and plot proportion of candidate and background genes 
## for GO-categories that have a FWER < 0.8
set.seed(123)
candi_gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59',
 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
bg_gene_ids = c('FGR', 'NPHP1', 'DRD2', 'ABCC10', 'PTBP2', 'JPH4',
 'SMARCC2', 'FN1', 'NODAL', 'CYP1A2', 'ACSS1', 'CDHR1', 'SLC25A36',
 'LEPR', 'PRPS2', 'TNFAIP3', 'NKX3-1', 'LPAR2', 'PGAM2', 'GAPDHS')
genes = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
names(genes) = c(candi_gene_ids, bg_gene_ids)
go_bg = go_enrich(genes, n_randsets=100, silent=TRUE)
plot_odds_ratio(go_bg, fwer_threshold=0.8)

## plot user-defined GO-IDs 
plot_odds_ratio(go_bg, go_ids=c('GO:0005634','GO:0004945',
 '0.05309471','GO:0008289','GO:0005737','GO:0071495'))

}

