\name{assignMutations}
\alias{assignMutations}
\title{Mutation Assignment}
\description{Assigns mutations to previously predicted subpopulations.}
\usage{assignMutations(dm, finalSPs, max_PM=6, cnvSPs=NULL, ploidy = 2, verbose = T)}
\arguments{
  \item{dm}{
       Matrix in which each row corresponds to a mutation. Has to contain at least the following column names:\cr 
	\bold{chr} - the chromosome on which each mutation is located; \cr
	\bold{startpos} - the genomic position of each mutation;\cr
	\bold{AF_Tumor} - the allele-frequency of each mutation;\cr
	\bold{PN_B} - the count of the B-allele in normal (non-tumor) cells (binary variable: 1 if the mutation is a germline variant, 0 if somatic).
}
  \item{finalSPs}{
       Matrix in which each row corresponds to a subpopulation, as calculated by \code{\link{clusterCellFrequencies}}.
}
  \item{max_PM}{
       Upper threshold for the number of amplicons per mutated cell. See also \code{\link{cellfrequency_pdf}}.
}
  \item{cnvSPs}{
       Matrix in which each row corresponds to a subpopulation, as calculated by \code{\link{clusterCellFrequencies}}. If not set, finalSPs will be used to assign CNVs as well as SNVs.
}
  \item{ploidy}{
        The background ploidy of the sequenced sample (default: 2). Changing the value of this parameter is not recommended. Dealing with cell lines or tumor biopsies of very   high (>=0.95) tumor purity is a necessary but not sufficient condition to change the value of this parameter.
}
\item{verbose}{
        Give a more verbose output.
}
}
\details{
Each mutated locus \eqn{l} is assigned to the subpopulation \eqn{C}, whose size \eqn{f_C} can best explain the allele frequency (AF) and copy number (CN) observed at \eqn{l}. Four alternative cell frequency probabilities, \eqn{P_x(f_C)}, are calculated for the SNV at locus \eqn{l}, with \eqn{x} denoting one of the four alternative evolutionary scenarios (see also \code{\link{cellfrequency_pdf}}).\cr
The SNV is assigned to subpopulation:\cr
\eqn{C:=argmax_C (P_s(f_C), P_p(f_C), P_c(f_C), P_i(f_C)) } (see \code{\link{cellfrequency_pdf}}).\cr

The mutated loci assigned to each subpopulation cluster represent the genetic profile of each predicted subpopulation.\cr
The assignment between subpopulation \eqn{C} and locus \eqn{l} only implies that the SNV at \eqn{l} has been first propagated during the clonal expansion that gave rise to \eqn{C}. So SNVs present in \eqn{C} may not be exclusive to \eqn{C} but may also be present in subpopulations smaller than \eqn{C}. Whether or not this is the case can sometimes be inferred from the phylogenetic structure of the subpopulation composition. See also \code{\link{buildPhylo}}.
  
}
\value{A list with two fields:
\item{dm}{
       The input matrix with seven additional columns:\cr
	\bold{SP} - subpopulation to which the point mutation has been assigned;\cr
	\bold{PM_B} - count of the B-allele at the mutated genomic locus, in the assigned subpopulation (SP).\cr
   	\bold{PM} - total count of all alleles, in the assigned subpopulation (SP).\cr
	\bold{SP_cnv} - if the point mutation lies within an amplified or deleted region: the subpopulation to which the copy number variation has been assigned.	This entry has the same value as SP if and only if: i) the SNV and the CNV were propagated during the same clonal expansion or ii) the SNV lies within a copy neutral region.\cr
	\bold{PM_B_cnv} - count of the B-allele, in the CNV harboring subpopulation (SP_cnv).\cr
	\bold{PM_cnv} - total count of all alleles, in the CNV harboring subpopulation (SP_cnv).\cr
	\bold{\%maxP} - confidence of the assigned SP/SP_cnv scenario.\cr
} 
\item{finalSPs}{
       The input matrix of subpopulations with column \bold{nMutations} updated according to the total number of mutations assigned to each subpopulation.
}
}

\references{Li, B. & Li, J. Z (2014). A general framework for analyzing tumor subclonality using SNP array and DNA sequencing data. Genome Biol.}

\author{Noemi Andor}


\seealso{
\code{\link{clusterCellFrequencies}}
}

