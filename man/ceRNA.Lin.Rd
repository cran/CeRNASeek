\name{ceRNA.Lin}
\alias{ceRNA.Lin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
identifying miRNA sponge interactions using ceRNA.Lin function 
}
\description{
identifying miRNA sponge interactions using ceRNA.Lin.In this function We
implement several popular linear methods ( HyperC, SC) to identify miRNA sponge.
interactions.
}
\usage{
ceRNA.Lin(miRtar, targetce = NULL, geneexp, miRexp, numMIR = 1, method = "pearson", 
numrandom = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{miRtar}{
A data frame representing the relationship between miRNA and target. The data frame contains the name of 
the miRNA and target regulatory relationship.Required option for method "pearson" and "partial correlation".
}
  \item{targetce}{
a character string (vector) specifying candidate target name to analyse (default (targetce = NULL)).
}
  \item{geneexp}{
An input target expression data frame, the columns are genes and
the rows are samples.The expression value may be  gene expression ,non-coding RNA expression  or 
expression values of circRNAs and so on. Required option for method pearson" and "partial correlation".
}
  \item{miRexp}{
An input miRNA expression data frame, the columns are miRNA and
the rows are samples. Required option for method pearson" and "partial correlation".
}
  \item{numMIR}{
a numeric vector that specify the minimum number of 2 gene-shared miRNAs.
}
  \item{method}{
a character string (default "pearson") indicating which statistical method to choose to calculate the ceRNA 
interaction relationship. One of "pearson" (default), or "partial correlation", can be abbreviated.
}
  \item{numrandom}{
The number of random. Required option for method "partial correlation",default
(numrandom = 100).
}
}
\details{
Note:All the arguments without default value must be assigned.
}
\value{
A list of identified miRNA sponge interactions containing following components:
\itemize{
\item \code{ceRNA} predicted triplets and related information,a 5 columns dataframe as following:
\itemize{
\item \code{targetce} represented target names,respectively.
\item \code{anotherce} names of modulators that another target(modulators) constitutes a ceRNA interaction relation.
\item \code{miRNAs} names of miRNA shared by two targets.
\item \code{miRNAs_num} number of miRNAs shared by two targets.
\item \code{pvalue} The p value of the identified ceRNA interaction relation.
}
\item \code{miR_l} Number of miRNAs interacting with each target in the input file.
}
}
\references{
Paci P , Colombo T , Farina L . Computational analysis identifies a sponge interaction network between
long non-coding RNAs and messenger RNAs in human breast cancer[J]. BMC Systems Biology, 2014, 8(1):83.
Zhang Y , Xu Y , Feng L , et al. Comprehensive characterization of lncRNA-mRNA related ceRNA network
across 12 major cancers[J]. Oncotarget, 2014, 7(39):64148-64167.
}
\examples{
##identifying miRNA sponge interactions
##Here we take the regulatory relationship between six genes and 71 miRNAs  and corresponding
##expression profilesas an example.
ceRNA.Lin(miRtar=dataset[["miRtar"]], targetce = NULL, geneexp=dataset[["geneexp"]],
          miRexp=dataset[["miRexp"]], numMIR = 1,method = "pearson", numrandom = 100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ceRNA.Lin}
