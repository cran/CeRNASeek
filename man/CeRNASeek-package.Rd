\name{CeRNASeek-package}
\alias{CeRNASeek-package}
\alias{CeRNASeek}
\docType{package}
\title{
CeRNASeek
}
\description{
identifying miRNA sponge interactions .
}
\details{
CeRNASeek considers five method to identify ceRNA crosstalk, including two types of global 
ceRNA regulation prediction methods(ratio based, we termed ratio and hypergeometric test
based, termed HyperT) and three types of context-specific prediction methods (hypergeometric 
test plus coexpression, termed HyperC, sensitivity correlation-based method (SC)and conditional 
mutual information (CMI)-based methods and cernia method).
   Ratio-based prediction:
this method ranked the candidate genes by the proportion of common miRNAs,For instance, we need to identify the 
ceRNA partners for gene i from all candidate gene sets S, and the ratio is calculated as the intersection of miRNAi 
and miRNAj divided by miRNAj Where miRNAi is the miRNA set that regulated gene i and miRNAj is the miRNA set that regulated gene j.

   Hypergeometric test-based prediction-HyperT:
it is usually used the hypergeometric to evaluate whether two genes were coregulated by miRNAs,This statistic 
test computed the significance of common miRNAs for each RNA pairs.

   Hypergeometric test combined with coexpression-based prediction-HyperC: 
To discover the active ceRNA-ceRNA regulatory pairs in a specific context, the commonly used method is to using the
coexpression principle to filter the ceRNA-ceRNA regulation identified based on the above two global methods,this
method integrated context-specific gene expression profile data sets. The Pearson correlation coefficient (R) of 
each candidate ceRNA regulatory pairs identified was calculated.

   SC-based prediction:
Another method introduced the expression profile data of shared miRNAs, and uses partial correlation to calculate ceRNA 
interaction pairs. Then, the Sensitivity correlation(SC) of miRNA-M, for the corresponding candidate ceRNA pair is 
calculated. 

   CMI-based methods:
CMI is widely used to identify the RNA-RNA correlations, given the value of miRNAs. Hermes is a widely used method, 
which predicts ceRNA interactions from expression profiles of candidate RNAs and their common miRNA regulators using CMI,
as in Sumazin et al.Firstly computed the significance of difference of I(miRNAi;T2|T1) and I(miRNAi;T2),where miRNAi represents the ith 
miRNA in the miRNA shared by the two target genes. Then Shuffled condition target's expression in 1000 times. The final 
significance based on Fisher’s method was calculated.For more information, please refer to the article by Sumazin P et al.

  cernia method
 cernia method is  implemented based on the following seven scores: 
1. The fraction of common miRNAs;
2. The density of the MREs for all shared miRNAs;
3. The distribution of MREs of the putative ceRNAs;
4. The relation between the overall number of MREs for a putative ceRNA, compared with the number of miRNAs that yield
these MREs;
5. The density of the hybridization energies related to MREs for
all shared miRNAs;
6. The DT-Hybrid recommendation scores; and
7. The pairwise Pearson correlation between putative ceRNA
expressions from selected tissue.
}
\references{
Xu J , Feng L , Han Z , et al. Extensive ceRNA–ceRNA interaction networks mediated by miRNAs regulate
development in multiple rhesus tissues[J]. Nucleic Acids Research, 2016:gkw587.

Sumazin P , Yang X , Chiu H S , et al. An Extensive MicroRNA-Mediated Network of RNA-RNA Interactions 
Regulates Established Oncogenic Pathways in Glioblastoma[J]. Cell, 2011, 147(2):0-381.

Paci P , Colombo T , Farina L . Computational analysis identifies a sponge interaction network between
long non-coding RNAs and messenger RNAs in human breast cancer[J]. BMC Systems Biology, 2014, 8(1):83.

Zhang Y , Xu Y , Feng L , et al. Comprehensive characterization of lncRNA-mRNA related ceRNA network
across 12 major cancers[J]. Oncotarget, 2014, 7(39):64148-64167.

Sardina D S , Alaimo S , Ferro A , et al. A novel computational method for inferring competing endogenous
interactions[J]. Briefings in Bioinformatics, 2016:bbw084.
}