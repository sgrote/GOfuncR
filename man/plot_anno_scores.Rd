\name{plot_anno_scores}
\alias{plot_anno_scores}
\title{
Plot distribution of scores of genes annotated to GO-categories
}
\description{
Uses the result of a GO-enrichment analysis performed with \code{\link{go_enrich}} and a vector of GO-IDs and plots for each of these GO-IDs the scores of the annotated genes. This refers to the scores that were provided as user-input in the \code{\link{go_enrich}} analysis.\cr
\code{plot_anno_scores} works with all four tests implemented in \code{\link{go_enrich}} (hypergeometric, wilcoxon rank-sum, binomial and 2x2 contingency table test), with test-specific output (see details).
}

\usage{
plot_anno_scores(res, go_ids, annotations = NULL)
}
\arguments{
    \item{res}{
        an object returned from \code{\link{go_enrich}} (list() of 3 elements).
    }
    \item{go_ids}{
        character() vector of GO-IDs, e.g. c('GO:0005737','GO:0071495').
        This specifies the GO-categories that are plotted.
    }
    \item{annotations}{
        optional data.frame() for custom annotations, with two columns: (1) gene-symbols and (2) GO-categories.
        This is needed if \code{\link{go_enrich}} was run with custom `annotations` to generate `res`
    }
}

\details{
The plot depends on the statistical test that was specified in the \code{\link{go_enrich}} call.

For the hypergeometric test pie charts show the amounts of candidate and background genes that are annotated to the GO-categories and the root nodes (candidate genes in the colour of the corresponding root node). The top panel shows the odds-ratio and 95\%-CI from fisher's exact test (two-sided) comparing the GO-categories with their root nodes. Note that \code{\link{go_enrich}} reports the the hypergeometric test for over- and under-representation of candidate genes wich correspond to the one-sided fisher's exact tests.

For the wilcoxon rank-sum test violin plots show the distribution of the scores of genes that are annotated to each GO-category and the root nodes. Horizontal lines in the left panel indicate the median of the scores that are annotated to the root nodes. The wilcoxon rank-sum test reported in the \code{\link{go_enrich}} result compares the scores annotated to a GO-category with the scores annotated to the corresponding root node.

For the binomial test pie charts show the amounts of A and B counts associated with each GO-category and root node, (A in the colour of the corresponding root node). The top-panel shows point estimates and the 95\%-CI of p(A) in the nodes, as well as horizontal lines that correspond to p(A) in the root nodes. The p-value in the returned object is based on the null hypothesis that p(A) in a node equals p(A) in the corresponding root node. Note that \code{\link{go_enrich}} reports that value for one-sided binomial tests.

For the 2x2 contingency table test pie charts show the proportions of A and B, as well as C and D counts associated with a GO-category. Root nodes are not shown, because this test is independent of the root category. The top panel shows the odds ratio and 95\%-CI from Fisher's exact test (two-sided) comparing A/B and C/D inside one node. Note that in \code{\link{go_enrich}}, if all four values are >=10, a chi-square test is performed instead of fisher's exact test.
}

\value{
For the hypergeometric, binomial and 2x2 contingency table test, a data.frame() with the statistics that are used in the plots.\cr
For the wilcoxon rank-sum test no statistical results are plotted, just the distribution of annotated scores. The returned element in this case is a list() with three data frames: annotations of genes to the GO-categories, annotations of genes to the root nodes and a table which contains for every GO-ID the corresponding root node.
}

\author{
Steffi Grote
}

\seealso{
\code{\link{go_enrich}}\cr
\code{\link{get_anno_genes}}\cr
\code{\link{get_names}}\cr
\code{\link{vioplot}}\cr
}

\examples{ 
#### see the package vignette for more examples

#### Note that argument 'n_randsets' is reduced 
#### to lower computational time in the example.

## Assign two random counts and to some random genes to create example input
set.seed(123)
high_A_genes = c('G6PD', 'GCK', 'GYS1', 'HK2', 'PYGL', 'SLC2A8',
    'UGP2', 'ZWINT', 'ENGASE')
low_A_genes = c('CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1',
    'GYG1', 'PAX2')
A_counts = c(sample(15:25, length(high_A_genes)),
    sample(5:15, length(low_A_genes)))
B_counts = c(sample(5:15, length(high_A_genes)),
    sample(15:25, length(low_A_genes)))
genes = data.frame(gene=c(high_A_genes, low_A_genes), A_counts, B_counts)

## perform enrichment analysis to find GO-categories with high fraction of A
go_binom = go_enrich(genes, test='binomial', n_randsets=20)

## plot sums of A and B counts associated with the top GO-categories
top_gos = head(go_binom[[1]]$node_id)
stats = plot_anno_scores(go_binom, go_ids=top_gos)

## look at the results of binomial test used for plotting
## (this is two-sided, go_enrich reports one-sided tests)
head(stats)

}


