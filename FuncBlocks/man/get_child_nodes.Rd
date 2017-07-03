\name{get_child_nodes}
\alias{get_child_nodes}
\title{
Get all child nodes of gene ontology categories
}
\description{
Returns all child nodes (sub-categories) of GO-categories given their GO-IDs, e.g. c('GO:0042254', 'GO:0000109'). The output also states the shortest distance to the child node. 
Note that a GO-ID itself is also considered as child with distance 0.
}
\usage{
get_child_nodes(go_ids)
}
\arguments{
  \item{go_ids}{
a vector of GO-IDs, e.g. c('GO:0051082', 'GO:0042254').
}
}

\value{
a dataframe with four columns: parent GO-ID, child GO-ID, child GO-name and distance.
}
\references{
[1] Ashburner, M. et al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25, 25–29.
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
child_nodes = get_child_nodes(c('GO:0090070', 'GO:0000109'))
child_nodes
}

