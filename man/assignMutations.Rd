\name{assignMutations}
\alias{assignMutations}
\title{Mutation Assignment}
\description{Assigns mutations to previously predicted subpopulations.}
\usage{assignMutations(dm, finalSPs, max_PM=6, peakselection='localsum')}
\arguments{
  \item{dm}{
       Matrix in which each row corrsponds to a mutation. Has to contain at least the following columnnames:\cr 
	\bold{chr} - the chromosome on which each mutation is located; \cr
	\bold{startpos} - the genomic position of each mutation;\cr
	\bold{AF_Tumor} - the allele-frequency of each mutation;\cr
	\bold{PN_B} - the ploidy of the B-allele in normal (non-tumor) cells (binary variable: 1 if the mutation is a germline variant, 0 if somatic).
}
  \item{finalSPs}{
       Matrix in which each row corresponds to a subpopulation, as computed by \code{\link{clusterCellFrequencies}.}
}
  \item{max_PM}{
       Upper threshold for the number of amplicons per mutated cel<F6>l (default: 6). See also \code{\link{cellfrequency_pdf}.}
}
  \item{peakselection}{
       Strategy used to select mutation specific cell-frequency probability peaks. Options: 'maximum','localsum'. 
}
}
\details{
Each mutated locus \eqn{l} is assigned to the subpopulation \eqn{C}, whose size \eqn{f_C} can best explain the allele frequency (AF) and copy number (CN) observed at \eqn{l}.Three alternative cell frequency probabilities, \eqn{P_x(f_C)}, are calculated for the SNV at locus \eqn{l}, with \eqn{x} denoting three alternative evolutionary scenarios (see also \code{\link{cellfrequency_pdf}}):\cr
	1. \eqn{x:=s} --> Separate fit of SNV and CNV. CNV does not influence ploidy of the SNV, either because CNV occurs before SNV or because SNV and CNV occur independently from each other (i.e. they are never co-propagated during the same clonal expansion)\cr
	2. \eqn{x:=p} --> Partial dependency of SNV ploidy on CNV. The SNV is propagated during the expansion of \eqn{C}. Subsequently, the CNV is propagated during a clonal expansion of a cell-member of \eqn{C}. \cr
	3. \eqn{x:=j} --> Joint fit of SNV and CNV, assuming they co-occur together in the same cell and are propagated during the exact same clonal expansion.\cr

If peakselection is set to 'maximum', then the SNV is assigned to subpopulation:\cr
\eqn{C:=argmax_C (P_s(f_C), P_p(f_C), P_j(f_C)) }.\cr
If peakselection is set to 'localsum', then the SNV is assigned to subpopulation:\cr
\eqn{C:=argmax_C (L(P_s, f_C), L(P_p, f_C), L(P_j, f_C)) },where:
\eqn{L(P,f_C):= \sum_{f \in peak(f_C)} P(f) } calculates the sum of all probabilities in individual peaks.\cr

The mutated loci assigned to each subpopulation cluster represent the genetic profile of each predicted subpopulation.\cr
The assignment between subpopulation \eqn{C} and locus \eqn{l} only implies that the SNV at \eqn{l} has been first propagated during the clonal expansion that gave rise to \eqn{C}. So SNVs present in \eqn{C} may not be exclusive to \eqn{C} but may also be present in subpopulations smaller than \eqn{C}. Whether or not this is the case can sometimes be inferred from the phylogenetic structure of the subpopulation composition. See also \code{\link{buildPhylo}}.
  
}
\value{A list with two fields:
\item{dm}{
       The input matrix with seven additional columns:\cr
	\bold{SP} - the subpopulation to which the point mutation has been asssigned;\cr
    \bold{PM_B} - the ploidy of the B-allele at the mutated genomic locus, in the assigned subpopulation (SP).\cr
   	\bold{PM} - the total ploidy of all alleles, in the assigned subpopulation (SP).\cr
	\bold{SP_cnv} - if the point mutation lies within an amplified or deleted region: the subpopulation to which the copy number variation has been asssigned.	This entry has the same value as SP if and only if: i) the SNV and the CNV were propagated during the same clonal expansion or ii) the SNV lies within a copy neutral region.\cr
	\bold{PM_cnv} - the total ploidy of all alleles, in the CNV harboring subpopulation (SP_cnv).\cr
	\bold{\%maxP} - confidence of the point mutation assignment to SP.\cr
	\bold{scenario} - the evolutionary scenario under which the subpopulation configurations for this genomic locus have been solved (see also parameter "snv_cnv_flag" in \code{\link{cellfrequency_pdf}}).
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

