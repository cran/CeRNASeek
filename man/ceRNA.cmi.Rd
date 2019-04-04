\name{ceRNA.cmi}
\alias{ceRNA.cmi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
identifying miRNA sponge interactions using ceRNA.cmi function. 
}
\description{
identifying miRNA sponge interactions using ceRNA.cmi.In this function We
implement CMI methods to identify miRNA sponge interactions.
}
\usage{
ceRNA.cmi(miRtar, targetce = NULL, geneexp, miRexp, numMIR = 1, num_perm = 100, 
cutoff = 0.05)
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
\item{geneexp}{
An input target expression data frame, the columns are genes and
the rows are samples.The expression value may be  gene expression ,non-coding RNA expression  or 
expression values of circRNAs and so on.
}
\item{miRexp}{
An input miRNA expression data frame, the columns are miRNA and
the rows are samples.
}
 \item{numMIR}{
a numeric vector that specify the minimum number of 2 gene-shared miRNAs
}
\item{num_perm}{
The number of random default(num_perm=100).
}
\item{cutoff}{
a numeric vector of the form method specifying threshold between
ceRNA interactions default(cutoff=0.05).
}
}
\details{
Note:All the arguments without default value must be assigned.The miRNA 
in the file of the target-miRNA regulatory relationship pair should also be 
present in the expression profile data file.
}
\value{
A list of identified miRNA sponge interactions containing following components:
\itemize{
\item \code{ceRNA_cmi} predicted triplets and related information,a 5 columns dataframe as following:
\itemize{
\item \code{targetce} represented target names,respectively.
\item \code{miRNA} names of miRNA in the triplet.
\item \code{anotherce} names of modulators that another target(modulators) constitutes a ceRNA interaction relation.
\item \code{cmi} Conditional mutual information(CMI) of triplets calculated using expression values.
\item \code{pvalue} The p value of the identified ceRNA interaction relation by CMI.
}
\item \code{ceRNA_comP} Number of miRNAs interacting with each target in the input file.
\itemize{
\item \code{targetce} represented target names,respectively.
\item \code{anotherce} names of predicted another target(modulators) constitutes a ceRNA interaction relation.
\item \code{commonMiR} names of miRNA shared by predicted ceRNA.
\item \code{xsq} Chi-square value of p-value \code{pvalue} in conditional mutual information .
\item \code{comP} The p value calculated by incomplete gamma function .
}
}
}
\references{
Sumazin P , Yang X , Chiu H S , et al. An Extensive MicroRNA-Mediated Network of RNA-RNA Interactions 
Regulates Established Oncogenic Pathways in Glioblastoma[J]. Cell, 2011, 147(2):0-381.
}
\examples{
##identifying miRNA sponge interactions.
##Here we take six candidate targets(modulators) and corresponding expression 
##data for example,Specify target(PTEN) to predict ceRNA interaction.
ceRNA.cmi(miRtar=dataset[["miRtar"]], targetce = "PTEN", geneexp=dataset[["geneexp"]],
          numMIR = 1, miRexp=dataset[["miRexp"]],num_perm=50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ceRNA.cmi}
