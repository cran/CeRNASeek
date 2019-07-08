\name{ceRNA.app}
\alias{ceRNA.app}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 identify the miRNAs that co-regulate the gene pairs of interest.
}
\description{
The ceRNA function is used to identify the miRNAs that co-regulate the gene pairs of interest.
The users can input the genes of interest. Otherwise, all the paired genes were calculated the
number of co-regulating miRNAs. 
}
\usage{
ceRNA.app(miRtar, targetce = NULL)
}
\arguments{
  \item{miRtar}{
A data frame representing the relationship between miRNA and target. The data frame contains the name of 
the miRNA and target regulatory relationship.
}
  \item{targetce}{
a character string (vector) specifying candidate target name to analyse (default (targetce = NULL)).
}
}
\details{
Note:All the arguments without default value must be assigned.
}
\value{
A list of identified  miRNAs that co-regulate the gene pairs containing following components:
\itemize{
\item \code{ceRNA}  identify the miRNAs that co-regulate the gene pairs of interest,List all possible 
ceRNA interactions,a 4 columns dataframe as following:
\itemize{
\item \code{targetce} represented target names,respectively.
\item \code{anotherce} names of modulators that another possible target(modulators) constitutes a ceRNA interaction relation.
\item \code{miRNAs} names of miRNA shared by two targets.
\item \code{miRNAs_num} number of miRNAs shared by two targets.
}
\item \code{miR_l} Number of miRNAs interacting with each target in the input file.
}
}

\examples{
##Here we take the regulatory relationship between six genes and 71 miRNAs.

ceRNA.app(dataset[["miRtar"]],targetce=NULL)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ceRNA.app}
