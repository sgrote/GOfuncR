\name{get_parent_nodes}
\alias{get_parent_nodes}
\title{
Get all parent nodes of gene ontology categories
}
\description{
Returns all parent nodes (higher level categories) of GO-categories given their GO-IDs, e.g. c('GO:0042254', 'GO:0000109'). The output also states the shortest distance to the parent node.
Note that a GO-ID itself is also considered as parent with distance 0.
}
\usage{
get_parent_nodes(go_ids)
}
\arguments{
  \item{go_ids}{
a vector of GO-IDs, e.g. c('GO:0051082', 'GO:0042254').
}
}

\value{
a dataframe with four columns: child GO-ID, parent GO-ID, parent GO-name and distance.
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
## get the parent nodes (higher level GO-categories) of three random GO-IDs
parents = get_parent_nodes(c('GO:0051082', 'GO:0042254', 'GO:0000109'))
parents
}

