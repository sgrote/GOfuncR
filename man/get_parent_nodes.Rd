\name{get_parent_nodes}
\alias{get_parent_nodes}
\title{
Get all parent nodes of gene ontology categories
}
\description{
Returns all parent nodes (higher level categories) of GO-categories given their GO-IDs, e.g. c('GO:0042254', 'GO:0000109').
The output also states the shortest distance to the parent node.
Note that a GO-ID itself is also considered as parent with distance 0.
}
\usage{
get_parent_nodes(go_ids, term_df = NULL, graph_path_df = NULL, godir = NULL)
}
\arguments{
    \item{go_ids}{
        a character() vector of GO-IDs, e.g. c('GO:0051082', 'GO:0042254').
    }
    \item{term_df}{
        optional data.frame() with an ontology 'term' table.
        Alternative to the default integrated GO-graph or \code{godir}.
        Also needs \code{graph_path_df}.
    }
    \item{graph_path_df}{
        optional data.frame() with an ontology 'graph_path' table.
        Alternative to the default integrated GO-graph or \code{godir}.
        Also needs \code{term_df}.
    }
    \item{godir}{
        optional character() specifying a directory that
        contains the ontology tables 'term.txt' and 'graph_path.txt'.
        Alternative to the default integrated GO-graph
        or \code{term_df} + \code{graph_path_df}.
    }
}

\details{
    By default the package's integrated GO-graph is used, but a custom ontology can be defined, too.
    For details on how to use a custom ontology with
    \code{term_df} + \code{graph_path_df} or \code{godir} please refer to the 
    package's vignette. The advantage of \code{term_df} + \code{graph_path_df}
    over \code{godir} is that the latter reads the files 'term.txt' and
    'graph_path.txt' from disk and therefore takes longer.
}

\value{
a data.frame() with four columns: child GO-ID (character()), parent GO-ID (character()), parent GO-name (character()) and distance (numeric()).
}
\references{
[1] Ashburner, M. et al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25, 25-29.
}
\author{
Steffi Grote
}

\seealso{
\code{\link{get_names}}\cr
\code{\link{get_child_nodes}}\cr
}

\examples{
## get the parent nodes (higher level GO-categories) of two random GO-IDs
parents = get_parent_nodes(c('GO:0051082', 'GO:0042254'))
parents
}

