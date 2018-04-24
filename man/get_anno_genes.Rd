\name{get_anno_genes}
\alias{get_anno_genes}
\title{
Get genes that are annotated to GO-categories
}
\description{
Given a vector of GO-IDs, e.g. c('GO:0072025','GO:0072221') the function returns all genes that are annotated to these GO-categories.
This includes genes that are annotated to any of the child nodes of a GO-category.
Besides the default 'Homo.sapiens', also other OrganismDb or OrgDb packages from Bioconductor, like 'Mus.musculus' (mouse) or 'org.Pt.eg.db' (chimp), can be used.
It is also possible to directly provide a dataframe with annotations, which is then searched for the input GO-categories and their child nodes.
By default the package's integrated GO-graph is used to find child nodes,
but a custom ontology can be defined, too.

}

\usage{
get_anno_genes(go_ids, database = 'Homo.sapiens', genes = NULL, annotations = NULL,
    term_df = NULL, graph_path_df = NULL, godir = NULL)
}

\arguments{
    \item{go_ids}{
        character() vector of GO-IDs, e.g. c('GO:0051082', 'GO:0042254').
    }
    \item{database}{
        optional character() defining an OrganismDb or OrgDb annotation package from Bioconductor, like 'Mus.musculus' (mouse) or 'org.Pt.eg.db' (chimp).
    }
    \item{genes}{
        optional character() vector of gene-symbols. If defined, only annotations of those genes are returned.
    }
    \item{annotations}{
        optional data.frame() with two columns: gene-symbols and GO-categories.
        Alternative to 'database'.
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
A data.frame() with two columns: GO-IDs (character()) and the annotated genes (character()). The output is ordered by GO-ID and gene-symbol.\cr 
}

\references{
[1] Ashburner, M. et al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25, 25-29.
}
\author{
Steffi Grote
}

\seealso{
\code{\link{get_anno_categories}}\cr
\code{\link{get_ids}}\cr
\code{\link{get_names}}\cr
\code{\link{get_child_nodes}}\cr
\code{\link{get_parent_nodes}}\cr
}

\examples{ 

## find all genes that are annotated to GO:0000109
## ("nucleotide-excision repair complex")
get_anno_genes(go_ids='GO:0000109')

## find out wich genes from a set of genes
## are annotated to some GO-categories
genes = c('AGTR1', 'ANO1', 'CALB1', 'GYG1', 'PAX2')
gos = c('GO:0001558', 'GO:0005536', 'GO:0072205', 'GO:0006821')
anno_genes = get_anno_genes(go_ids=gos, genes=genes)
# add the names and domains of the GO-categories
cbind(anno_genes ,get_names(anno_genes$go_id)[,2:3])

## find all annotations to GO-categories containing 'serotonin receptor'
sero_ids = get_ids('serotonin receptor')
sero_anno = get_anno_genes(go_ids=sero_ids$go_id)
# merge with names of GO-categories
head(merge(sero_ids, sero_anno))
}

