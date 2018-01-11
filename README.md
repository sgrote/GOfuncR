### GOfuncR: Gene Ontology Enrichment Using FUNC 

_GOfuncR_ performs a gene ontology enrichment analysis based on the ontology enrichment software [_FUNC_](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1800870/). It provides the standard candidate vs. background enrichment analysis using the hypergeometric test, as well as three additional tests: (i) the Wilcoxon rank-sum test that is used when genes are ranked, (ii) a binomial test that can be used when genes are associated with two counts, e.g. amino acid changes since a common ancestor in two different species, and (iii) a Chi-square or Fisher's exact test that is used in cases when genes are associated with four counts, e.g. non-synonymous or synonymous variants that are fixed between or variable within species.  
To correct for multiple testing and interdependency of the tests, family-wise error rates (FWER) are computed based on random permutations of the gene-associated variables. _GOfuncR_ also provides tools for exploring the ontology graph and the annotations, and options to take gene-length or spatial clustering of genes into account during testing.  
GO-annotations and gene-coordinates are obtained from _OrganismDb_ packages ([_Homo.sapiens_](https://www.bioconductor.org/packages/release/data/annotation/html/Homo.sapiens.html) by default) or _OrgDb_ and _TxDb_ packages. The gene ontology graph (obtained from [geneontology](http://archive.geneontology.org/latest-termdb/), last modified on 12-Oct-2017), is integrated in the package. 


#### Installation

##### from Bioconductor

_GOfuncR_ is currently available only from the ['devel' version of Bioconductor](https://bioconductor.org/packages/3.7/bioc/html/GOfuncR.html). 
It will be in the stable 'release' version when Bionductor 3.8 is released.


##### from GitHub

Optional vignette (tutorial) generation needs _pandoc_ installed. 
There may also be other dependencies that need to be installed manually.

+ Installation with vignette


```r
## install packages needed for vignette generation in R
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocStyle","Homo.sapiens"))
```

```bash
## install GOfuncR from the console
git clone https://github.com/sgrote/GOfuncR
R CMD build GOfuncR
R CMD INSTALL GOfuncR_0.99.11.tar.gz
```


+ Installation without vignette

The vignette can also be found [here](https://bioconductor.org/packages/3.7/bioc/vignettes/GOfuncR/inst/doc/GOfuncR.html) and does not need to be installed.

```bash
## install GOfuncR from the console
git clone https://github.com/sgrote/GOfuncR
R CMD build --no-build-vignettes GOfuncR
R CMD INSTALL GOfuncR_0.99.11.tar.gz
```

or

```r
## install from GitHub in R:
install.packages("devtools")
library(devtools)
install_github("sgrote/GOfuncR")
```

#### Usage  

See the [vignette](https://bioconductor.org/packages/3.7/bioc/vignettes/GOfuncR/inst/doc/GOfuncR.html) on Bioconductor for an introduction.

If _GOfuncR_ was installed with vignette, the vignette can also be opened in R:
```r
## open tutorial in R
browseVignettes("GOfuncR")
```

Also see the man-pages for single functions, e.g.
```r
library(GOfuncR)
?go_enrich
?get_anno_categories

```
