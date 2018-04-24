\name{get_anno_categories}
\alias{get_anno_categories}
\title{
Get all associated ontology categories for the input genes
}
\description{
Returns all associated GO-categories given a vector of gene-symbols, e.g. c('SPAG5', 'BTC'). 
Besides the default 'Homo.sapiens', also other OrganismDb or OrgDb packages from Bioconductor, like 'Mus.musculus' (mouse) or 'org.Pt.eg.db' (chimp), can be used.
It also possible to directly provide a dataframe with annotations, which is then searched for the input genes and filtered for GO-categories that are present in the ontology.
}
\usage{
get_anno_categories(genes,  database = 'Homo.sapiens', annotations = NULL, termdf = NULL, godir = NULL, silent = FALSE)
}
\arguments{
    \item{genes}{
        a character() vector of gene-symbols, e.g. c('SPAG5', 'BTC').
    }
    \item{database}{
        optional character() defining an OrganismDb or OrgDb annotation package from Bioconductor, like 'Mus.musculus' (mouse) or 'org.Pt.eg.db' (chimp).
    }
    \item{annotations}{
        optional data.frame() with two columns: gene-symbols and GO-categories.
        Alternative to 'database'.
    }
    \item{termdf}{
        optional data.frame() with at an ontology 'term' table.
        Alternative to the default integrated GO-graph or \code{godir}.
    }
    \item{godir}{
        optional directory (string) that contains 'term.txt' table.
        Alternative to the default integrated GO-graph or \code{termdf}.
    }
    \item{silent}{
        logical. If TRUE all output to the screen except for warnings and errors is suppressed.
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
Also note that GO-categories which are not represented or obsolete in the internal GO-graph of GOfuncR or the custom ontology provided through \code{termf} or \code{godir} are removed to be consistent with the annotations used in \code{\link{go_enrich}}.
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

