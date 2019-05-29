\name{refine}
\alias{refine}
\title{
Refine results given a FWER threshold 
}
\description{
Given a FWER threshold, this function refines the results from \code{go_enrich()} like described in the elim algorithm of [1].\cr
This algorithm removes genes from significant child-categories and then checks whether a category is still significant.\cr
This way significant results are restricted to more specific categories.

}
\usage{
    refine(res, fwer = 0.05, fwer_col = 7, annotations = NULL)
}
\arguments{
    \item{res}{
        list() returned from \code{go_enrich()} 
    }
    \item{fwer}{
        numeric() FWER-threshold. Categories with a FWER < \code{fwer} will be labeled significant.
    }
    \item{fwer_col}{
        numeric() or character() specifying the column of \code{go_enrich()[[1]]} that is to be filtered.
        E.g. 6 for under-representation or 7 over-representation of candidate genes in the hypergeometric test.
    }
    \item{annotations}{
        optional data.frame() with custom annotations. Only needed if \code{go_enrich()} was run with custom annotations in the first place.
    }
}

\details{
    For each domain a p-value is found by interpolation, that corresponds to the input FWER threshold.
    Since GO-domains are independent graphs, the same FWER will correspond to different p-values, e.g. in 'molecular_function' and 'biological_process'.
}

\value{
a data.frame() with p-values after refinement for categories that were significant in \code{go_enrich()[[1]]} given the FWER-threshold.
}
\references{
[1] Alexa, A. et al. (2006). Improved scoring of functional groups from gene expression data by decorrelating GO graph structure. Bioinformatics 22, 1600-1607.

}
\author{
Steffi Grote
}



\seealso{

\code{\link{go_enrich}}\cr

}
\examples{

## perform enrichment analysis for some genes
gene_ids = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1',
    'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
input_hyper = data.frame(gene_ids, is_candidate=1)
res_hyper = go_enrich(input_hyper, n_randset=100, silent=TRUE)
head(res_hyper[[1]])
## perform refinement for categories with FWER < 0.1
refined = refine(res_hyper, fwer=0.1)
refined

}



