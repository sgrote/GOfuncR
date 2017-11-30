\name{get_anno_categories}
\alias{get_anno_categories}
\title{
Get all associated ontology categories for the input genes
}
\description{
Returns all associated GO-categories given a vector of gene-symbols, e.g. c('SPAG5', 'BTC'). 
Besides the default 'Homo.sapiens', also other OrganismDb or OrgDb packages from Bioconductor, like 'Mus.musculus' (mouse) or 'org.Pt.eg.db' (chimp), can be used.
}
\usage{
get_anno_categories(genes,  database = 'Homo.sapiens')
}
\arguments{
    \item{genes}{
        a character() vector of gene-symbols, e.g. c('SPAG5', 'BTC').
    }
    \item{database}{
        optional character() defining an OrganismDb or OrgDb annotation package from Bioconductor, like 'Mus.musculus' (mouse) or 'org.Pt.eg.db' (chimp).
    }
}

\value{
a data.frame() with two columns: gene (character()) and GO-ID (character()).
}
\references{
[1] Ashburner, M. et al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25, 25-29.
}
\author{
Steffi Grote
}

\note{
This gives only direct annotations of genes to GO-categories. By definition genes are also indirectly annotated to all parent nodes of those categories. Use \code{\link{get_parent_nodes}} to get the higher level categories of the directly annotated GO-categories.\cr
Also note that GO-categories which are not represented or obsolete in the internal GO-graph of GOfuncR are removed to be consistent with the annotations used in \code{\link{go_enrich}}.
}

\seealso{
\code{\link{get_anno_genes}}\cr
\code{\link{get_parent_nodes}}\cr
\code{\link{get_names}}\cr
}

\examples{

## get the GO-annotations for two random genes
anno1 = get_anno_categories(c('BTC', 'SPAG5'))
head(anno1)
# add the names and domains of the GO-categories
named_anno1 = merge(anno1, get_names(anno1$go_id))
head(named_anno1)

}

