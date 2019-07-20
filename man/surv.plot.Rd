\name{surv.plot}
\alias{surv.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Draw survival curve
}
\description{
It is used to draw the survival curve of ternary relationship pairs 
and to support the survival prognosis of training sets and test sets.
}
\usage{
surv.plot(ceRNA, exp.sur, train = NULL, test = NULL, index)
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
Survival curve of user-selected ceRNA ternary pairs.
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
surv.plot(ceRNA=dataset[["Pre.ceRNA"]],exp.sur=dataset[["exp.sur"]],train=NULL,test=NULL,index=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ceRNA.surv}