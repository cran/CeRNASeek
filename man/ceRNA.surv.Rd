\name{ceRNA.surv}
\alias{ceRNA.surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
survival analysis of ceRNA ternary relationship pairs
}
\description{
It is used to predict the survival curve of ternary relationship pairs 
and to support the survival prognosis of training sets and test sets.
}
\usage{
ceRNA.surv(ceRNA, exp.sur, train = NULL, test = NULL, index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ceRNA}{
a dataframe that the ceRNA relationship is the data, and the prediction result data obtained 
according to the ceRNA prediction algorithm.Such as ceRNA.cmi prediction result file.
}
  \item{exp.sur}{
 dataframe specifying expression and survival information.Its rownames are sample names.Its colnames
are names in triplets and survival information (see example data in details).
}
  \item{train}{
a character string vector specifying train sample names.
}
  \item{test}{
a character string vector specifying test sample names.
}
  \item{index}{
a numeric vector (default 1) reprsenting rowindex of triplets analyed.
}
}

\value{
A list of identified miRNA sponge interactions containing following components:
\itemize{
\item \code{targetce} represented target names,respectively.
\item \code{anotherce} names of modulators that another target(modulators) constitutes a ceRNA interaction relation.
\item \code{coef_targetce} the \emph{coxph} coefficient of targetce;
\item \code{p_targetce} the \emph{coxph} significance of targetce;
\item \code{coef_anotherce}  the \emph{coxph} coefficient of anotherce;
\item \code{p_anotherce} the \emph{coxph} significance of anotherce;
\item \code{N_low} The genes were ranked according to expression, and the number of samples 
                    expressed in the bottom 25%.
\item \code{N_high} The genes were ranked according to expression, and the number of samples 
                    expressed in the top 25%.
\item \code{HR} Risk score of ceRNA ternary pair.
\item \code{p} Survival significance of ceRNA ternary pairs.
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
ceRNA.surv(ceRNA=dataset[["Pre.ceRNA"]],exp.sur=dataset[["exp.sur"]],train=NULL,test=NULL,index=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ceRNA.surv}