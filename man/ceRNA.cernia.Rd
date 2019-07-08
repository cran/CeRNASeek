\name{ceRNA.cernia}
\alias{ceRNA.cernia}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
identifying miRNA sponge interactions using ceRNA.cernia function.
}
\description{
identifying miRNA sponge interactions using ceRNA.cernia.In this function We
implement cernia score to identify miRNA sponge interactions.
The fraction of common miRNAs;
1. The fraction of common miRNAs;
2. The density of the MREs for all shared miRNAs;
3. The distribution of MREs of the putative ceRNAs;
4. The relation between the overall number of MREs for a putative ceRNA, 
   compared with the number of miRNAs that yield these MREs;
5. The density of the hybridization energies related to MREs for
   all shared miRNAs;
6. The DT-Hybrid recommendation scores; 
7. The pairwise Pearson correlation between putative ceRNA
   expressions.

}
\usage{
ceRNA.cernia(miRtar, targetce = NULL, geneexp, miRexp, mres, numMIR = 3, cor_cutoff = 0,
             s_cutoff = 0.5)
}
\arguments{
  \item{miRtar}{
A data frame representing the relationship between miRNA and target. The data frame contains the name of 
the miRNA and target regulatory relationship.
}
  \item{targetce}{
a character string (vector) specifying candidate target name to analyse (default (targetce = NULL)).
}
  \item{geneexp}{
An input target expression data frame, the columns are genes and
the rows are samples.The expression value may be  gene expression ,non-coding RNA expression  or 
expression values of circRNAs and so on.
}
  \item{miRexp}{
An input miRNA expression data frame, the columns are miRNA and
the rows are samples.
}
  \item{mres}{
 miRNA response elements (mres) data frame, each row
contains five elements: mirna, target, energy, gap_l, gap_r. 
}
  \item{numMIR}{
a numeric vector that specify the minimum number of 2 gene-shared miRNAs.
}
  \item{cor_cutoff}{
a numeric vector of the form method specifying threshold between
ceRNA interactions correlation,default(cor_cutoff=0).
}
  \item{s_cutoff}{
the threshold of seven comprehensive scores. 
}
}
\details{
Note:All the arguments without default value must be assigned.The miRNA 
in the file of the target-miRNA regulatory relationship pair should also be 
present in the expression profile data file.
Internal functions (parMM, graphWeights, recommendation, dtHybrid) of cernia method are from
the website: https://github.com/dsardina/cernia Copyright 2016 Rosalba Giugno Licensed under
the Apache License, Version 2.0 (the 'License')
}
\value{
A list of identified miRNA sponge interactions containing following components:
\itemize{
\item \code{targetce} represented target names,respectively.
\item \code{anotherce} names of modulators that another target(modulators) constitutes a ceRNA interaction relation.
\item \code{miRNAs_num} names of miRNA in the triplet.
\item \code{Score 1} The fraction of common miRNAs;
\item \code{Score 2} The density of the MREs for all shared miRNAs;
\item \code{Score 3} The distribution of MREs of the putative ceRNAs;
\item \code{Score 4} The relation between the overall number of MREs for a putative ceRNA, 
                     compared with the number of miRNAs that yield these MREs;
\item \code{Score 5} The density of the hybridization energies related to MREs for
                     all shared miRNAs;
\item \code{Score 6} The DT-Hybrid recommendation scores; 
\item \code{Score 7} The pairwise Pearson correlation between putative ceRNA
                     expressions.
\item \code{Normalized score} combine the sum of 7 scores and standardize by z-score.



}
}
\references{
 Sardina D S , Alaimo S , Ferro A , et al. A novel computational method for inferring 
 competing endogenous interactions[J]. Briefings in Bioinformatics, 2016:bbw084.
 
}

\examples{
##identifying miRNA sponge interactions.
##Here we take six candidate targets(modulators) and corresponding expression 
##data for example,Specify target(PTEN) to predict ceRNA interaction.
ceRNA.cernia(miRtar=dataset[["miRtar"]], targetce = "PTEN", geneexp=dataset[["geneexp"]],
          numMIR = 1, miRexp=dataset[["miRexp"]], mres=dataset[["mres"]], 
		  cor_cutoff = 0.2, s_cutoff = 0.5) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ceRNA.cernia}
