\name{get_anno_categories}
\alias{get_anno_categories}
\title{
Get all associated ontology categories for the input genes
}
\description{
Returns all associated GO-categories given a vector of gene-symbols, e.g. c('SPAG5', 'BTC').
}
\usage{
get_anno_categories(genes, ref_genome = 'grch37')
}
\arguments{
  \item{go_ids}{
a vector of gene-symbols, e.g. c('SPAG5', 'BTC').
}
  \item{ref_genome}{
One of 'grch37', 'grch38' or 'grcm38'. This specifies the reference genome used for GO-annotations.
}
}

\value{
a dataframe with four columns: gene, GO-ID, GO-name and GO-domain.
}
\references{
[1] Ashburner, M. et al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25, 25–29.
}
\author{
Steffi Grote
}

\seealso{
\code{\link{get_anno_genes}}\cr
\code{\link{get_parent_nodes}}\cr
}

\examples{
## get the GO-annotations for two random genes
get_anno_categories(c('BTC', 'SPAG5'))

## get the GO-annotations for two random mouse genes
get_anno_categories(c('Mus81', 'Papola'), ref_genome='grcm38')
}

