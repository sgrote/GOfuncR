\name{plot_scores}
\alias{plot_scores}
\title{
Plot distribution of scores for genes in enriched or user-defined GO-categories
}
\description{
Uses the result of a wilcoxon rank-sum test performed with \code{\link{go_enrich}} and extracts the GO-categories that are significantly (given a FWER-threshold) enriched in genes with high scores. Alternatively GO-IDs can be defined explicitly.\cr
\code{plot_scores} then creates a violin plot for each GO-category showing the distribution of the scores of annotated genes. These are the scores that were provided as input for the wilcoxon rank-sum test. Above each violin plot the number of annotated genes is stated. The colors of the violin plots indicate the GO-domain.\cr
The plot also shows the distribution of scores in the three root nodes (domains) of the gene ontology.

}

\usage{
plot_scores(res, fwer_threshold = 0.05, go_ids = NULL)
}
\arguments{
  \item{res}{
an object returned from \code{\link{go_enrich}} with the option \code{test='wilcoxon'} (list of 3 elements).
}
  \item{fwer_threshold}{
numeric defining the FWER-threshold. \code{plot_scores} extracts all GO-categories from \code{res} that have a FWER < \code{fwer_threshold} and plots the distribution of the scores of annotated genes for those GO-categories. Only used when \code{go_ids} are not specified.
}
  \item{go_ids}{
optional vector of GO-IDs, e.g. c('GO:0005737','GO:0071495'). If defined, this specifies the GO-categories directly and the \code{fwer_threshold} is not used.
}
}

\details{
The plotted scores correspond to the input for the wilcoxon rank-sum test, i.e. \code{res[[2]]}.
}

\value{
A list with two data frames, 'anno_score' and 'anno_score_root' which contain the genes and scores annotated to each of the plotted GO-categories and the root nodes (invisible).
}

\author{
Steffi Grote
}


\seealso{
\code{\link{go_enrich}}\cr
\code{\link{get_anno_genes}}\cr
\code{\link{get_names}}\cr
\code{\link{vioplot}}\cr
\code{\link{plot_odds_ratio}}\cr
}

\examples{ 
#### Note that argument 'n_randsets' is reduced 
#### to lower computational time in the examples.

## perform a wilcoxon rank-sum test for some random genes
## and plot distribution of annotated scores 
## for GO-categories that have a FWER < 0.21
set.seed(123)
high_score_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8', 'UGP2', 'ZWINT', 'ENGASE')
low_score_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = c(sample(20:30, length(high_score_genes)), sample(5:15, length(low_score_genes)))
names(genes) = c(high_score_genes, low_score_genes)
go_willi = go_enrich(genes, test='wilcoxon', n_randsets=100)
plot_scores(go_willi, fwer_threshold=0.21)


## plot user-defined GO-IDs 
plot_scores(go_willi, go_ids=c('GO:0005634','GO:0004945',
 'GO:0008289','GO:0005737','GO:0071495'))


}

