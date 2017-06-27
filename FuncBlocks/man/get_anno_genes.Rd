\name{get_anno_genes}
\alias{get_anno_genes}
\title{
Get genes that are annotated to enriched or user-defined GO-categories
}
\description{
Uses an object returned from \code{\link{go_enrich}} as input and returns the GO-categories that are significantly (given a FWER-threshold) enriched in candidate genes (or genes with high scores if \code{test='wilcoxon'}, together with the genes that are annotated to those GO-categories.
Alternatively, also user-defined GO-categories can be used as input.
}

\usage{
get_anno_genes(res, fwer_threshold = 0.05, background = FALSE, 
  go_ids = NULL, ref_genome = NULL, genes = NULL)
}
\arguments{
  \item{res}{
an object returned from \code{\link{go_enrich}} (list of 3 elements). If not defined, \code{go_ids} have to be specified. 
}
  \item{fwer_threshold}{
numeric defining the FWER-threshold. Given \code{res} as input, \code{get_anno_genes} extracts all GO-categories from \code{res} that have a FWER < \code{fwer_threshold} and adds the (candidate) genes that are annotated to those GO-categories.  
}
  \item{background}{
logical indicating whether background genes should be included. Only used when \code{res} is defined and contains the results from a hypergeometric test (which is the default in \code{\link{go_enrich}}). 
}
  \item{go_ids}{
vector of GO-IDs, e.g. c('GO:0051082', 'GO:0042254'). If \code{res} is not defined, \code{go_ids} specifies the GO-categories for which annotated genes will be returned.
}
  \item{ref_genome}{
One of 'grch37', 'grch38' or 'grcm38'. If \code{res} is not defined, this specifies the reference genome used for GO-annotations.
}
  \item{genes}{
optional vector of gene-symbols. If defined, only annotations of those genes are returned. Only used when \code{res} is not defined.
}
}

\value{
A data frame with the GO-IDs and the annotated genes for significantly enriched or user-defined GO-categories. If \code{res} was given as input, two additional columns are added: the FWER and the score which was used for the genes in the \code{\link{go_enrich}} input (1/0 for candidate and background genes for the hypergeometric test or scores for the wilcoxon rank sum test).\cr
Given \code{res} as input, the output is ordered by the FWER, else it is ordered by GO-ID.\cr 
}

\references{
[1] Ashburner, M. et al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25, 25–29.
}
\author{
Steffi Grote
}


\seealso{
\code{\link{go_enrich}}\cr
\code{\link{get_names}}\cr
\code{\link{get_child_nodes}}\cr
\code{\link{get_parent_nodes}}\cr
}

\examples{ 
#### Note that argument 'n_randsets' is reduced 
#### to lower computational time in the examples.
 
## perform an enrichment analysis for some candidate genes
## and get candidate genes annotated to GO-categories that have a FWER < 0.05
set.seed(123)
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2',
 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
genes = rep(1, length(gene_ids))
names(genes) = gene_ids
go_res = go_enrich(genes, n_randset=100, silent=TRUE)
get_anno_genes(go_res, fwer_threshold=0.05)

## find all genes that are annotated to GO:0000109
## ("nucleotide-excision repair complex")
get_anno_genes(go_ids="GO:0000109", ref_genome="grch38")

}

