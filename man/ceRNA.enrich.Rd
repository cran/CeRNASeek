\name{ceRNA.enrich}
\alias{ceRNA.enrich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Enrich the biological functions of interest. 
}
\description{
Downstream analysis function of ceRNA,This function can be used to identify 
biological functions of interest, and users can enrich the functions of 
interest by ceRNA.enrich.
}
\usage{
ceRNA.enrich(data, GOterms, background, threshold = 2, correction = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a dataframe that the ceRNA relationship is the data, and the prediction result data obtained 
according to the ceRNA prediction algorithm.Such as ceRNA.cmi prediction result file.
}
  \item{GOterms}{
a list that goterm of interest and all gene sets in the term.
}
  \item{background}{
a vector containing a gene set in which GOterm annotated genes must be.Its id style must be 
consistent with the id format in \emph{GOterms}.
}
  \item{threshold}{
a numeric (default 2) representing min number of intersection between a modulator's targets and a GOterms genes.
}
  \item{correction}{
correction method (default "BH") in one of \code{\link{p.adjust.methods}}.
}
}
\details{
Note:All the arguments without default value must be assigned.
}
\value{
A list of identified miRNA sponge interactions containing following components:
\itemize{
\item \code{target} represented target names,respectively.
\item \code{GOterm} the GOterm name.
\item \code{target_num} names of represented target in the triplet.
\item \code{GOtermnum} the gene number of a GOterm;
\item \code{term_tar} the number of intersected factor between a GOterm genes and a modulator targets;
\item \code{P_value} the p value of the significance enrichment;
\item \code{fdr} corrected P_value by the assigned method;

}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
ceRNA.enrich(data=dataset[["Pre.ceRNA"]],GOterms=dataset[["GOterms"]],
             background=dataset[["background"]],threshold=1,correction="BH")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ceRNA.enrich}