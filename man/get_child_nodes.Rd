\name{get_child_nodes}
\alias{get_child_nodes}
\title{
Get all child nodes of gene ontology categories
}
\description{
Returns all child nodes (sub-categories) of GO-categories given their GO-IDs, e.g. c('GO:0042254', 'GO:0000109').
The output also states the shortest distance to the child node. 
By default the package's integrated GO-graph is used, but a custom ontology can be defined, too.
Note that a GO-ID itself is also considered as child with distance 0.


}
\usage{
get_child_nodes(go_ids, term_df = NULL, graph_path_df = NULL, godir = NULL)
}
\arguments{
    \item{go_ids}{
        a character() vector of GO-IDs, e.g. c('GO:0051082', 'GO:0042254').
    }
    \item{term_df}{
        optional data.frame() with at an ontology 'term' table.
        Alternative to the default integrated GO-graph or \code{godir}.
        Also needs \code{graph_path_df}.
    }
    \item{graph_path_df}{
        optional data.frame() with at an ontology 'graph_path' table.
        Alternative to the default integrated GO-graph or \code{godir}.
        Also needs \code{term_df}.
    }
    \item{godir}{
        optional directory (string) that contains 'term.txt'
        and 'graph_path.txt' table.
        Alternative to the default integrated GO-graph
        or \code{term_df} + \code{graph_path_df}.
    }
}

\details{
    For details on how to use a custom ontology with
    \code{term_df} + \code{graph_path_df} or \code{godir} please refer to the 
    package's vignette. The advantage of \code{term_df} + \code{graph_path_df}
    over \code{godir} is that the latter reads the files 'term.txt' and
    'graph_path.txt' from disk and therefore takes longer.
}

\value{
a data.frame() with four columns: parent GO-ID (character()), child GO-ID (character()), child GO-name (character()) and distance (numeric()).
}
\references{
[1] Ashburner, M. et al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25, 25-29.
}
\author{
Steffi Grote
}

\seealso{
\code{\link{get_names}}\cr
\code{\link{get_parent_nodes}}\cr
}

\examples{
## get the child nodes (sub-categories) of two random GO-IDs
child_nodes = get_child_nodes(c('GO:0090070', 'GO:0000112'))
child_nodes
}

