\name{ceRNA.basic}
\alias{ceRNA.basic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying miRNA sponge interactions using ceRNA.basic function 
}
\description{
This function predicts ceRNA interactions by ratio or hypergeometric test.
}
\usage{
ceRNA.basic(miRtar, targetce = NULL, method = "ratio", numMIR = 1,
 cutoff = ifelse(method == "ratio", 1/3, 0.05))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{miRtar}{
A data frame representing the relationship between miRNA and target. The data frame contains the name of 
the miRNA and target regulatory relationship.
}
  \item{targetce}{
a character string (vector) specifying candidate target name to analyse (default (targetce = NULL)).
}
  \item{method}{
a character string (default "ratio") indicating which statistical method to choose to calculate the ceRNA 
interaction relationship. One of "ratio" (default), or "hypergeometric", can be abbreviated.
}
  \item{numMIR}{
a numeric vector that specify the minimum number of 2 gene-shared miRNAs
}
  \item{cutoff}{
a numeric vector of the  method("ratio","hypergeometric") specifying threshold between ceRNA interactions .(default (1/3)).
}
}
\details{
Note:All the arguments without default value must be assigned.
}
\value{
A list of identified miRNA sponge interactions containing following components:
\itemize{
\item \code{cesig} predicted  significant triplets and related information,a 5 columns dataframe as following:
\itemize{
\item \code{targetce} represented target names,respectively.
\item \code{anotherce} names of modulators that another target(modulators) constitutes a ceRNA interaction relation.
\item \code{miRNAs} names of miRNA shared by two targets.
\item \code{miRNAs_num} number of miRNAs shared by two targets.
\item \code{ratio/pvalue} The ratio/pvalue(optional for method("ratio","hypergeometric")) of the identified ceRNA interaction relation.
}
\item \code{cenotsig} predicted  not significant triplets and related information,a 5 columns dataframe as following:
\itemize{
\item \code{targetce} same as \code{targetce} in \code{cesig}
\item \code{anotherce} same as \code{anotherce} in \code{cesig}
\item \code{miRNAs} same as \code{miRNAs} in \code{cesig}
\item \code{miRNAs_num} same as \code{miRNAs_num} in \code{cesig}
\item \code{pvalue} same as \code{pvalue} in \code{cesig}
}
}
}
\references{
Xu J , Feng L , Han Z , et al. Extensive ceRNA–ceRNA interaction networks mediated by miRNAs regulate
development in multiple rhesus tissues[J]. Nucleic Acids Research, 2016:gkw587.
}
\examples{
##identifying miRNA sponge interactions
##Here we take six candidate targets(modulators) for example
ceRNA.basic(miRtar=dataset[["miRtar"]],targetce=NULL,method="ratio",numMIR = 1,cutoff =1/3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ceRNA.basic}
