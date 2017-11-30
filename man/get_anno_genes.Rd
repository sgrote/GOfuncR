\name{get_anno_genes}
\alias{get_anno_genes}
\title{
Get genes that are annotated to GO-categories
}
\description{
Given a vector of GO-IDs, e.g. c('GO:0072025','GO:0072221') the function returns all genes that are annotated to these GO-categories.
This includes genes that are annotated to any of the child nodes of a GO-category.
Besides the default 'Homo.sapiens', also other OrganismDb or OrgDb packages from Bioconductor, like 'Mus.musculus' (mouse) or 'org.Pt.eg.db' (chimp), can be used.
It is also possible to only consider a user-defined set of genes.
}

\usage{
get_anno_genes(go_ids, database = 'Homo.sapiens', genes = NULL)
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

