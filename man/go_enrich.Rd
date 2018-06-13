\name{go_enrich}
\alias{go_enrich}
\title{
Test gene sets for enrichment in GO-categories
}
\description{
Tests GO-categories for enrichment of user defined gene sets, using either the hypergeometric (default), 
Wilcoxon rank-sum, binomial or 2x2 contingency table test.

}

\usage{
go_enrich(genes, test = 'hyper', n_randsets = 1000, organismDb = 'Homo.sapiens', gene_len = FALSE,
    regions = FALSE, circ_chrom = FALSE, silent = FALSE, domains = NULL, orgDb = NULL,
    txDb = NULL, annotations = NULL, gene_coords = NULL,  godir = NULL)
}

\arguments{

\item{genes}{
a data.frame() with gene-symbols (character()) in the first column and test-dependent additional numeric() columns:\cr
If \code{test='hyper'} (default) a second column with 1 for candidate genes and 0 for background genes. 
If no background genes are defined, all remaining genes from the internal dataset are used as background.
All candidate genes are implicitly part of the background gene set.\cr
If \code{test='wilcoxon'} a second column with the score that is associated with each gene.\cr
If \code{test='binomial'} two additional columns with two gene-associated integers.\cr
If \code{test='contingency'} four additional columns with four gene-associated integers.\cr
For \code{test='hyper'} with \code{regions=TRUE} the first column describes chromosomal regions ('chr:start-stop', e.g. '9:0-39200000').\cr
Note that each gene has to be unique in the data.frame; e.g. for \code{test='wilcoxon'} a gene cannot have two different scores assigned.
}

\item{test}{
    character(). 'hyper' (default) for the hypergeometric test, 'wilcoxon' for the Wilcoxon rank-sum test, 'binomial' for the binomial test and 'contingency' for the 2x2-contingency table test (fisher's exact test or chi-square).
}

\item{n_randsets}{
    integer defining the number of random sets for computing the FWER.
}  

\item{organismDb}{
    character(). Annotation package for GO-annotations and gene coordinates.
    Besides the default 'Homo.sapiens' also 'Mus.musculus' and 'Rattus.norvegicus' are currently available on Bioconductor.
}

\item{gene_len}{
    logical. Correct for gene length. If \code{TRUE} the probability of a background gene to be chosen as a candidate gene in a random set is dependent on the gene length.
    If \code{FALSE} genes are chosen randomly with equal probability each. Only for \code{test='hyper'}.
}

\item{regions}{
    logical. If \code{TRUE} chromosomal regions are analyzed, and genes are automatically extracted from the regions defined in \code{genes[,1]} as e.g. '9:0-39200000'. Note that this option requires the input of background regions.
}

\item{circ_chrom}{
    logical. When \code{genes} defines chromosomal regions, \code{circ_chrom=TRUE} uses background regions from the same chromosome and allows randomly chosen blocks to overlap multiple background regions. Only if \code{test='hyper'}.
}

\item{silent}{
    logical. If \code{TRUE} all output to the screen except for warnings and errors is suppressed.
}

\item{domains}{
    optional character() vector containing one or more of the three GO-domains 'cellular_component', 'biological_process' and 'molecular_function'.
    If defined, the analysis will be reduced to those domains which saves time.\ cr
    If a custom ontology is specified in \code{godir} it might have a different domains.
}

\item{orgDb}{
    optional character() naming an OrgDb annotation package from Bioconductor. If \code{orgDb} is set, \code{organismDb} is not used. OrgDb annotations are available for a wider range of species, e.g. 'org.Pt.eg.db' for chimp and 'org.Gg.eg.db' for chicken. Note that options \code{gene_len} and \code{regions} also need a \code{txDb} for the gene coordinates (which are integrated in OrganismDb).
}

\item{txDb}{
    optional character() naming an TxDb annotation package from Bioconductor (e.g. 'TxDb.Ptroglodytes.UCSC.panTro4.refGene') for the gene coordinates.
    Only needed when \code{orgDb} is specified, and options \code{gene_len} or \code{regions} are set.
    Note that \code{orgDb} is needed whenever \code{txDb} is defined, even when custom \code{annotations} are provided, because the \code{orgDb} package is used for Entrez-ID to gene-symbol conversions.
}

\item{annotations}{
    optional data.frame() for custom annotations, with two character() columns:
    (1) gene-symbols and (2) GO-categories.
    Note that options \code{gene_len} and \code{regions} also need an \code{organismDb} or \code{txDb} + \code{orgDb}.
}

\item{gene_coords}{
    optional data.frame() for custom gene coordinates, with four columns:
    gene-symbols (character), chromosome (character), start (integer), end (integer).
    When \code{gene_len=TRUE} or \code{regions=TRUE} these custom gene coordinates are used instead of the ones obtained from \code{organismDb} or \code{txDb}.
}

\item{godir}{
    optional character() specifying a directory () that contains a custom GO-graph
    (files 'term.txt', 'term2term.txt' and 'graph_path.txt').
    Alternative to the default integrated GO-graph.\cr
    For details please refer to the package's vignette.
}


}

\details{
Please also refer to the package's vignette.\cr
GO-annotations are taken from a Bioconductor annotation package (OrganismDb package 'Homo.sapiens' by default), but also other 'OrganismDb' or 'OrgDb' packages can be used. It is also possible to provide custom annotations as a data.frame().

The ontology graph is integrated, but a custom version can be defined as well with parameter 'godir'.
As long as the ontology tables are in the right format (see link to description in vignette),
any ontology can be used in GOfuncR, it is not restricted to the gene ontology. 

The statistical analysis is based on the ontology enrichment software FUNC [2].
\code{go_enrich} offers four different statistical tests: (1) the hypergeometric test for a candidate and a background gene set; (2) the Wilcoxon rank-sum test for genes that are ranked by scores (e.g. p-value for differential expression); (3) the binomial test for genes that have two associated counts (e.g. amino-acid changes on the human and the chimp lineage); and (4) a 2x2-contingency table test for genes that have four associated counts (e.g. for a McDonald-Kreitman test).

To account for multiple testing family-wise error rates are computed using randomsets.
Besides naming candidate genes explicitly, for the hypergeometric test it is also possible to provide entire genomic regions as input.
The enrichment analysis is then performed for all genes located in or overlapping these regions and the multiple testing correction accounts for the spatial clustering of genes.
}

\value{
A list with components
\item{results}{
a data.frame() with the FWERs from the enrichment analyses per ontology category, ordered by 'FWER_overrep', 'raw_p_overrep', -'FWER_underrep', -'raw_p_underrep', 'ontology' and 'node_id', or the corresponding columns if another test then the hypergeometric test is used.
This table contains the combined results for all three ontology domains. Note that GO-categories without any annotations of candidate or background genes are not listed.}
\item{genes}{
the input \code{genes} data.frame(), excluding those genes for which no GO-annotations are available and which therefore were not included in the enrichment analysis. If \code{gene_len=TRUE}, also genes without gene coordinates are excluded.}
\item{databases}{
a data.frame() with the used annotation packages and GO-graph.}
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
\code{\link{plot_anno_scores}}\cr
\code{\link{get_names}}\cr
\code{\link{get_ids}}\cr
}

\examples{

#### Note that argument 'n_randsets' is reduced 
#### to lower computational time in the following examples. 
#### Using the default value is recommended.

#### Perform a GO-enrichment analysis for some human genes 
#### with a defined background set
# create input dataframe that defines the candidate and backround genes
candi_gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 
    'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
bg_gene_ids = c('FGR', 'NPHP1', 'DRD2', 'ABCC10', 'PTBP2', 'JPH4', 'SMARCC2',
    'FN1', 'NODAL', 'CYP1A2', 'ACSS1', 'CDHR1', 'SLC25A36', 'LEPR', 'PRPS2',
    'TNFAIP3', 'NKX3-1', 'LPAR2', 'PGAM2')
is_candidate = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
genes = data.frame(gene_ids=c(candi_gene_ids, bg_gene_ids), is_candidate)
genes

# run enrichment analysis
go_res = go_enrich(genes, n_randset=100)

# go_enrich returns a list with 3 elements:
# 1) results from the anlysis
# (ordered by FWER for overrepresentation of candidate genes)
head(go_res[[1]])
# see the top GOs from every GO-domain
by(go_res[[1]], go_res[[1]][,'ontology'], head)
# 2) all valid input genes
go_res[[2]]
# 3) annotation databases used
go_res[[3]]

#### see the package's vignette for more examples
}

\keyword{ htest }
