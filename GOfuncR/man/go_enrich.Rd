\name{go_enrich}
\alias{go_enrich}
\title{
Test human or mouse genes for enrichment in GO-categories
}
\description{
Tests all GO-categories from the three domains 'biological process', 'cellular component' and 'molecular function' for enrichment of user defined gene sets [1]. The statistical analysis is based on the ontology enrichment software FUNC [2].\cr
\code{go_enrich} offers four different statistical tests: (1) the hypergeometric test for a candidate and a background gene set; (2) the wilcoxon rank-sum test for genes that are ranked by scores (e.g. p-value for differential expression); (3) the binomial test for genes that have two associated counts (e.g. amino-acid changes on the human and the chimp lineage); and (4) a 2x2-contingency table test for genes that have four associated counts (e.g. for a McDonald-Kreitman test). 
Besides naming candidate genes explicitly, for the hypergeometric test it is also possible to provide entire genomic regions as input. The enrichment analysis is then performed for all genes located in or overlapping these regions and the multiple testing correction accounts for the spatial clustering of genes. Note that this option requires the input of background regions.\cr
The gene ontology, the annotations and the gene-coordinates are integrated in the package and updated regularly.
}
\usage{
go_enrich(genes, test = 'hyper', n_randsets = 1000, gene_len = FALSE, 
  circ_chrom = FALSE, ref_genome = 'grch37', silent = FALSE, domains = NULL)
}
\arguments{
  \item{genes}{
A dataframe with gene-symbols in the first column and test-dependent additional columns:\cr
If \code{test='hyper'} (default) a second column with 1 for candidate genes and 0 for background
genes. If no background genes are defined, all remaining genes from the internal dataset are used as background.\cr
If \code{test='wilcoxon'} a second column with the score that is associated with each gene.\cr
If \code{test='binomial'} two additional columns with two gene-associated integers.\cr
If \code{test='contingency'} four additional columns with four gene-associated integers.\cr
For \code{test='hyper'} the first column may also describe chromosomal regions ('chr:start-stop').
}
  \item{test}{
'hyper' (default) for the hypergeometric test, 'wilcoxon' for the Wilcoxon rank test, 'binomial' for the binomial test and 'contingency' for the 2x2-contingency table test (fisher's exact test or chi-square).
}
  \item{n_randsets}{
integer defining the number of random sets created to compute the FWER.
}  
  \item{gene_len}{
logical. If \code{test='hyper'} the probability of a background gene to be chosen as a candidate gene in a random set is dependent on the gene length. 
}
  \item{circ_chrom}{
logical. When \code{genes} defines chromosomal regions, \code{circ_chrom=TRUE} uses background regions from the same chromosome and allows randomly chosen blocks to overlap multiple background regions. Only if \code{test='hyper'}.
}
  \item{ref_genome}{
'grch37' (default), 'grch38' or 'grcm38'. Defines the reference genome used for GO-annotations and gene coordinates.
}
  \item{silent}{
logical. If TRUE all output to the screen except for warnings and errors is suppressed.
}
  \item{domains}{
optional vector containing one or more of the three GO-domains 'cellular_component', 'biological_process' and 'molecular_function'. If defined, the analysis will be reduced to those domains which saves time.
}

}

\details{
For details please refer to \code{vignette("GOfuncR", package="GOfuncR")}.
}

\value{
A list with components
\item{results}{
a dataframe with the FWERs from the enrichment analyses per ontology category, ordered by 'FWER_overrep', 'raw_p_overrep', -'FWER_underrep', -'raw_p_underrep', 'ontology' and 'node_id', or the corresponding columns if another test then the hypegeometric test is used. This table contains the combined results for all three ontology domains. Note that GO-categories without any annotations of candidate or background genes are not listed.}
\item{genes}{
the input \code{genes} dataframe, excluding those genes for which no GO-annotations are available and which therefore were not included in the enrichment analysis. If \code{gene_len=TRUE}, also genes without gene coordinates are excluded.}
\item{ref_genome}{
the reference genome used.}
}

\references{
[1] Ashburner, M. et al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25: 25-29. \doi{10.1038/75556}\cr
[2] Pruefer, K. et al. (2007). FUNC: A package for detecting significant associations between gene
sets and ontological. BMC Bioinformatics 8: 41. \doi{10.1186/1471-2105-8-41}
}
\author{
Steffi Grote
}

\seealso{
\code{\link{get_parent_nodes}}\cr
\code{\link{get_child_nodes}}\cr
\code{\link{get_anno_categories}}\cr
\code{\link{get_anno_genes}}\cr
\code{\link{plot_anno_scores_ratio}}\cr
\code{\link{get_names}}\cr
\code{\link{get_ids}}\cr
}

\examples{
#### see the package vignette for more examples

#### Note that argument 'n_randsets' is reduced 
#### to lower computational time in the following examples. 
#### Using the default value is recommended.

#### Perform a GO-enrichment analysis for some random human genes
# create input dataframe 

# TODO: change to data-frame, maybe add another test

gene_ids = c('NCAPG', 'QUATSCH1', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 
 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
genes
# run enrichment analysis
go_res = go_enrich(genes, n_randset=100)
# go_enrich returns a list with 3 elements:
# 1) results from the anlysis (ordered by FWER for overrepresentation of candidate genes)
head(go_res[[1]])
# see the top GOs from every GO-domain
by(go_res[[1]], go_res[[1]][,'ontology'], head)
# 2) all valid input genes
go_res[[2]]
# 3) reference genome used (default=grch37)
go_res[[3]]

}

\keyword{ htest }
