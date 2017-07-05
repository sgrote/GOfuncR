\name{go_enrich}
\alias{go_enrich}
\title{
Test human or mouse genes for enrichment in GO-categories
}
\description{
Tests all GO-categories from the three domains 'biological process', 'cellular component' and 'molecular function' for enrichment of user defined candidate genes compared to a set of background genes [1]. The statistical analysis is based on the ontology enrichment software FUNC [2].\cr
Besides naming candidate genes explicitly, it is also possible to provide entire genomic regions as input. The enrichment analysis is then performed for all genes located in or overlapping these regions and the multiple testing correction accounts for the spatial clustering of genes. Note that this option requires the input of background regions.\cr
The gene ontology, the annotations and the gene-coordinates are integrated in the package and updated regularly.
}
\usage{
go_enrich(genes, test = 'hyper', n_randsets = 1000, gene_len = FALSE, 
  circ_chrom = FALSE, ref_genome = 'grch37', silent = FALSE, domains = NULL)
}
\arguments{
  \item{genes}{
If \code{test='wilcoxon'} a numeric vector of scores.
If \code{test='hyper'} (default) a binary vector with 1 for candidate genes and 0 for background
genes. If no background genes are defined, all remaining genes from the internal dataset are used as background.
The names of the vector are the corresponding gene-symbols. 
For \code{test='hyper'} the names of the vector can also describe chromosomal regions ('chr:start-stop').
}
  \item{test}{
'hyper' (default) for the hypergeometric test or 'wilcoxon' for the Wilcoxon rank test.
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
%\details{
%For details please refer to \code{vignette("ABAEnrichment", package="ABAEnrichment")}.
%}
\value{
A list with components
\item{results}{
a dataframe with the FWERs from the enrichment analyses per ontology category, ordered by 'FWER_overrep', 'raw_p_overrep', -'FWER_underrep', -'raw_p_underrep', 'ontology' and 'node_id'. This table contains the combined results for all three ontology domains. Note that GO-categories without any annotations of candidate or background genes are not listed.}
\item{genes}{
a vector of the requested genes, excluding those genes for which no GO-annotations are available and which therefore were not included in the enrichment analysis. If \code{gene_len=TRUE}, also genes without gene coordinates are excluded.}
\item{ref_genome}{
the reference genome used.}   
}

\references{
[1] Ashburner, M. et al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25: 25–29. \doi{10.1038/75556}\cr
[2] Pruefer, K. et al. (2007). FUNC: A package for detecting significant associations between gene
sets and ontological. BMC Bioinformatics 8: 41. \doi{10.1186/1471-2105-8-41}
}
\author{
Steffi Grote
}

\seealso{
\code{\link{get_parent_nodes}}\cr
\code{\link{get_child_nodes}}\cr
\code{\link{get_anno_genes}}\cr
\code{\link{plot_odds_ratio}}\cr
\code{\link{get_names}}\cr
\code{\link{get_ids}}\cr
}

\examples{

#### Note that argument 'n_randsets' is reduced 
#### to lower computational time in the following examples. 
#### Using the default value is recommended.

#### Perform a GO-enrichment analysis for some human genes
# create input vector
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


#### Perform a GO-enrichment analysis for some random mouse genes
#### and account for gene-length
#### and consider only root nodes 'molecular function' and 'biological process'
# create input vector
gene_ids = c('Arsi', 'Mapk4', 'Papola', 'Tfrc', 'Bak1', 'Fopnl', 'Mus81', 'Opa3', 'Npcd')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
genes
# run enrichment analysis
go_mouse = go_enrich(genes, n_randsets=100, ref_genome='grcm38',
 gene_len=TRUE, domains=c('molecular_function','biological_process'))
# look at top results
head(go_mouse[[1]])


#### Perform a GO-enrichment analysis for some random genomic regions
# create input vector with names of the for 'chr:start:stop'
genes = c(1, rep(0,6))
names(genes) = c('8:81000000-83000000', '7:1300000-56800000', '7:74900000-148700000',
 '8:7400000-44300000', '8:47600000-146300000', '9:0-39200000', '9:69700000-140200000')
# one candidate and 6 background regions:
genes
# run enrichment analysis for genes located in these regions
go_region = go_enrich(genes, n_randsets=100)
# look at top results
head(go_region[[1]])
# genes located in the candidate region
go_region[[2]][go_region[[2]]==1]


# For more examples please refer to 
# /home/steffi_grote/R_packages/FuncBlocks_usage_examples.R.

}

\keyword{ htest }
