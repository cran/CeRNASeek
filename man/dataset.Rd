\name{dataset}
\alias{dataset}
\docType{data}
\title{
Data for Examples
}
\description{
This object contains data for examples.
}
\usage{data(dataset)}
\format{
A list with 4 variables:
  \describe{
	 \item{miRtar}{
	  MiRNA–target regulation,assembled genome-wide miRNA-gene regulation by TargetScan,It includes 44 rows standing for miRNA-gene regulation pairs in glioblastoma and 2 columns (the first is dataframe header).
	}\item{geneexp}{
	  a dataframe representing target regulations expression ,all of which are gene expression value,in glioblastoma.It includes 6 rows standing for 6 gene targets and 541 columns standing for 541 samples.Its rownames are gene symbols.
    }\item{miRexp}{
      a dataframe representing miRNA expression profile in glioblastoma.It includes 57 rows standing for 57 miRNA and 541 columns standing for 541 samples.Its rownames are miRNA symbols.
    }\item{Pre.ceRNA}{
      a data frame representing the ceRNA identified in the pre-experiment using the method provided by the software to draw a network map.
    }\item{mres}{
      a data frame representing the MRE.
    }\item{GOterms}{
      a list representing that goterm of interest and all gene sets in the term.
    }\item{background}{
      a character representing a gene set in which GOterm annotated genes must be.
    }\item{exp.sur}{
      a data frame representing specifying expression and survival information.Its rownames are sample names.Its colnames
      are names in triplets and survival information
    }\item{train}{
      a data frame representing specifying train sample names
    }\item{test}{
      a data frame representing specifying test sample names.
    }
  }
}
\details{
All epression data is from a study about glioblastoma.The miRNA-gene regulation information is from TargetScan.
}
\source{
Network, Atlas T C G. et al.
}
\references{
Network, Atlas T C G . Corrigendum: Comprehensive genomic characterization defines human glioblastoma genes and core pathways[J]. Nature, 2013, 494(7438):506-506.
}
\examples{
data(dataset)
## maybe str(dataset) ; ...
}
\keyword{datasets}

