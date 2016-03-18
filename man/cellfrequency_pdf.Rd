\name{cellfrequency_pdf}
\alias{cellfrequency_pdf}
\title{Computes the probability distribution of cellular frequencies for a single mutation.}
\description{Calculates \eqn{P} - the probability density distribution of cellular frequencies for one single point mutation or CNV. For each cell-frequency \eqn{f}, the value of \eqn{P(f)} reflects the probability that the mutation is present in a fraction \eqn{f} of cells.
}
\usage{cellfrequency_pdf(af, cnv, pnb, freq, max_PM=6, snv_cnv_flag=3, SP_cnv = NA, PM_cnv = NA)}
\arguments{
  \item{af}{
	  The allelic frequency at which the point muation has been observed.
}
  \item{cnv}{
	  The average copy number of the locus in which the mutation is embedded.
}
  \item{pnb}{
	  The ploidy of the B-allele in normal cells (binary variable: 1 if the mutation is a germline variant, 0 if somatic). B-alleles with normal cell ploidy>1 are not modeled.
}
  \item{freq}{
	  Array of cellular frequencies at which the probabilities will be calculated.
}
  \item{max_PM}{
	  Upper threshold for the number of amplicons per mutated cell (default: 6). \eqn{max\_PM } is the maximum number of amplicons above which solutions are rejected in the cell-frequency estimation step described below, i.e. \eqn{PM} <= \eqn{max\_PM. } The choice of \eqn{max\_PM }should depend on genomic depth of coverage and on the fraction of the genome sequenced: the higher the quality and abundance of data, the higher \eqn{max\_PM.}
}
  \item{snv_cnv_flag}{
 	  Flag indicating the evolutionary scenario under which frequency should be estimated: 1 - cellular frequency of SNV only; 2 - cellular frequency of CNV only (parameters AF and pnb are ignored); 3 - cellular frequency of SNV and CNV simultaneously, under the assumption that they are identical (co-occurrence assumption of SNV and CNV in same subpopulation); 4 - cellular frequency of SNV, under the assumption that it occured in the ancestor of the subpopulation with the CNV (this is the only scenario under which the SNV can have distinct ploidies in distinct subpopulations). Default: 3.
}
  \item{SP_cnv}{
	  Size of the subpopulation that harbors a copy number variation (CNV) at this locus. This variable is only relevant if an CNV and an SNV have overlapping genomic location, yet have been propagated during distinct clonal expansions (snv_cnv_flag=4).
}
  \item{PM_cnv}{
          Total ploidy in subpopulation which harbors a copy number variation (CNV) at this locus. This variable is only relevant if an CNV and an SNV have overlapping genomic location, yet have been propagated during distinct clonal expansions (snv_cnv_flag=4).
}
}
\details{
	We consider two types of molecular mechanisms that convert a locus into its mutated state: copy number variation (CNV) inducing events and single nucleotide variation (SNV) inducing events. We assume that a normal state is defined by a total ploidy of two and B allele ploidy below two, whereas a mutated state has an increased fraction of B alleles. The conditions defining these states for each locus \eqn{l} are as follows:\cr
	\bold{i}) \eqn{PM^B, PN^B, PM, PN \in N}; \bold{ii})  \eqn{PM^B \geq 1; PN^B \leq 1; PN = 2}; \bold{iii}) \eqn{\frac{PM^B}{PM} \geq \frac{PN^B}{PN} }. \cr

\eqn{PM^B} and \eqn{PN^B}  denote the ploidy of the B allele in each cell type: mutated cells and normal cells, respectively. The value of \eqn{PN^B} is one if \eqn{l} has a germline variant, zero otherwise. \eqn{PM, PN}  are the total ploidy of mutated cells and normal cells.  \eqn{PM} is required to be between one and \eqn{max\_PM}, that is, we exclude solutions for which the maximum number of amplicons per cell exceeds the user defined constant \eqn{max\_PM}. \cr
The function returns the probability distribution, \eqn{P(f)}, that the mutation at locus \eqn{l} is present in a fraction \eqn{f} of cells, where \eqn{f \in [0,1]}. Four alternative cell frequency probability distributions, \eqn{P(f)}, can be obtained for each allele-frequency + copy number pair (AF, CN).\cr

\bold{1.} \eqn{P_s(f_{cnv})} separately modeling the size \eqn{f_{cnv}} of the subpopulation propagating an CNV: \eqn{PM * f_{cnv} + PN *(1-f_{cnv}) = CN} \cr
2. \eqn{P_s(f_{snv})} and \eqn{P_p(f_{snv})} modeling the size \eqn{f_{snv}} of the subpopulation propagating an SNV: \cr
    \bold{2a)} \eqn{P_s(f_{snv})}: \eqn{PM^B * f_{snv} + PN^B *(1-f_{snv}) = AF*CN}, where \eqn{PM^B \leq max(2, PM)}; \cr
Here \eqn{f_{snv}} is calcualted separately of \eqn{f_{cnv}}, under the assumption that i) SNV and CNV occur independently from each other (i.e. they are never co-propagated during the same clonal expansion) or ii) SNV occured in a denscendant of the subpopulation with the CNV.\cr
    \bold{2b)}  \eqn{P_p(f_{snv})}: \eqn{PM^B * (f_{snv}-f_{cnv}) + pm^B * f_{cnv} + PN^B *(1-f_{snv}) = AF*CN}, where \eqn{pm^B \neq PM^B} and \eqn{pm^B \neq 2}.\cr
Here \eqn{f_{snv}} is calcualted partially dependent on \eqn{f_{cnv}}, under the assumption that the SNV occured in an ancestor of the subpopulation with the CNV. \cr
\bold{3.} \eqn{P_j(f)} jointly modeling the size \eqn{f} of the subpopulation propagating both the SNV and the CNV simulataneously: enforcing both equations, 1) and 2a), with additional constraints: \eqn{f:=f_{snv}=f_{cnv}} and \eqn{PM^B \leq PM} \cr

In 1) and 2) the SNV is present in a subpopulation different of the CNV harboring subpopulation. In 3) both the SNV and an CNV at \eqn{l} were propagated during the same clonal expansion.\cr 
}
\value{List with four components:
	  \item{p}{
	         The probability that the point mutation/CNV is present in a fraction \eqn{f} of cells, for each input frequency \eqn{f} in parameter \eqn{freq}.
}         
	  \item{bestF}{
                 The cellular frequency that best explains the observed allele frequency and/or copy number.
}
	   \item{fit}{
                Matrix with each row containing one alternative solution, (PM, PM_B, f), as well as an assesment of how well the solution fits above equations (Column "dev"). 
}
	   \item{errors}{
                Errors encountered during the density estimation step. 
}
}
\references{Noemi Andor, Julie Harness, Sabine Mueller, Hans Werner Mewes and Claudia Petritsch. (2013) ExPANdS: Expanding Ploidy and Allele Frequency on Nested Subpopulations. Bioinformatics.}
\author{Noemi Andor}
\examples{
freq=seq(0.1,1.0,by=0.01);
cfd=cellfrequency_pdf(af=0.26,cnv=1.95,pnb=0,freq=freq, max_PM=6)
plot(freq,cfd$p,type="l",xlab="f",ylab="P(f)");
}
