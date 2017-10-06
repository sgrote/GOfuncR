\name{get_ids}
\alias{get_ids}
\title{
Get the ID of a GO-category given its name
}
\description{
Returns GO-categories given (part of) their name. Matching is not case-sensitive.
}
\usage{
get_ids(go_name)
}
\arguments{
    \item{go_name}{
        (partial) name of a GO-category 
    }
}


\value{
a data frame with the full names of the GO-categories that contain \code{go_name}; together with the GO-domain ('cellular_component', 'biological_process' or 'molecular_function') and the GO-category IDs.
}
\references{
[1] Ashburner, M. et al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25, 25-29.

}
\author{
Steffi Grote
}

\note{
This is just a \code{grep(..., ignore.case=TRUE)} on the node names of the ontology. \cr
More sophisticated searches, e.g. with regular expressions, could be performed on the table returned by \code{get_ids('')} which lists all non-obsolete GO-categories.
}

\seealso{
\code{\link{get_names}}\cr
\code{\link{get_parent_nodes}}\cr
\code{\link{get_child_nodes}}\cr

}
\examples{
## get GO-IDs of categories that contain 'gabaergic' in their names
get_ids('gabaergic')

## get GO-IDs of categories that contain 'blood-brain barrier' in their names
get_ids('blood-brain barrier')

## get all valid GO categories
all_nodes = get_ids('')
head(all_nodes)
}



