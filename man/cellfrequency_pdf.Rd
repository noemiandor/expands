\name{cellfrequency_pdf}
\alias{cellfrequency_pdf}
\title{Computes the probability distribution of cellular frequencies for a single mutation.}
\description{Calculates \eqn{P} - the probability density distribution of cellular frequencies for one single point mutation or CNV. For each cell-frequency \eqn{f}, the value of \eqn{P(f)} reflects the probability that the mutation is present in a fraction \eqn{f} of cells.
}
\usage{cellfrequency_pdf(af, cnv, pnb, freq, max_PM=6, ploidy = 2)}
\arguments{
  \item{af}{
	  The allelic frequency at which the point mutation has been observed.
}
  \item{cnv}{
	  The average copy number of the locus in which the mutation is embedded.
}
  \item{pnb}{
	  The count of the B-allele in normal cells (binary variable: 1 if the mutation is a germline variant, 0 if somatic). B-alleles that have >1 copy in normal cells are not modeled.
}
  \item{freq}{
	  Vector of cellular frequencies at which the probabilities will be calculated.
}
  \item{max_PM}{
	  Upper threshold for the number of amplicons per mutated cell. \eqn{max\_PM } is the maximum number of amplicons above which solutions are rejected in the cell-frequency estimation step described below, i.e. \eqn{PM} <= \eqn{max\_PM. } The choice of \eqn{max\_PM }should depend on genomic depth of coverage and on the fraction of the genome sequenced: the higher the quality and abundance of data, the higher \eqn{max\_PM.}
}
  \item{ploidy}{
        The background ploidy of the sequenced sample (default: 2). Changing the value of this parameter is not recommended. Dealing with cell lines or tumor biopsies of very  high (>=0.95) tumor purity is a necessary but not sufficient condition to change the value of this parameter.
}
}
\details{
	We consider two types of molecular mechanisms that convert a locus into its mutated state: copy number variation (CNV) inducing events and single nucleotide variation (SNV) inducing events. We assume that a normal state is defined by a total allele count of two and B allele count below two, whereas a mutated state has an increased fraction of B alleles. The conditions defining these states for each locus \eqn{l} are as follows:\cr
	\bold{i}) \eqn{PM_B, PN_B, PM, PN \in N}; \bold{ii})  \eqn{PM_B \geq 1; PN_B \leq 1; PN = 2}; \bold{iii}) \eqn{\frac{PM_B}{PM} \geq \frac{PN_B}{PN} }. \cr

\eqn{PM_B} and \eqn{PN_B}  denote the count of the B allele in each cell type: mutated cells and normal cells, respectively. The value of \eqn{PN_B} is one if \eqn{l} has a germline variant, zero otherwise. \eqn{PM, PN}  are the total allele count of mutated cells and normal cells.  \eqn{PM} is required to be between one and \eqn{max\_PM}, that is, we exclude solutions for which the maximum number of amplicons per cell exceeds the user defined value of \eqn{max\_PM}. \cr
The function returns the probability distribution, \eqn{P(f)}, that the mutation at locus \eqn{l} is present in a fraction \eqn{f} of cells, where \eqn{f \in [0,1]}. \cr\cr

Four alternative cell frequency probability distribution scenarios, \eqn{P(f)}, can be obtained for each allele-frequency + copy number pair (AF, CN). For each scenario, model starts with a germline population that will be the root of all other modeled subpopulations. First subpopulation (\eqn{f_{cnv}}) modeled to evolve from the germline population is always the one carrying a CNV:\cr
\eqn{pm * f_{cnv} + PN *(1-f_{cnv}) = CN}, where \eqn{pm} is the total allele count of \eqn{f_{cnv}}.\cr
A subsequent subpopulation (\eqn{f_{snv}}) is always defined by an SNV and is modeled in relation to \eqn{f_{cnv}}, either as:\cr
\bold{1.} \eqn{P_s(f)} - its sibling: \eqn{PM_B * f_{snv} + PN_B *(1-f_{snv}) = AF*CN}, where \eqn{f_{snv}+f_{cnv}<=1; PM_B<=2}.\cr
\bold{2.} \eqn{P_p(f)} - its parent: \eqn{PM_B * (f_{snv}-f_{cnv}) + pm_B * f_{cnv} + PN_B *(1-f_{snv}) = AF*CN}, where \eqn{f_{snv}>f_{cnv}; PM_B<=2} and  \eqn{pm_B} is the B-allele count of \eqn{f_{cnv}}.\cr
\bold{3.} \eqn{P_c(f)} - its child: \eqn{PM_B * f_{snv} + PN_B *(1-f_{snv}) = AF*CN}, where \eqn{f_{snv}<f_{cnv}; PM_B<=pm}.\cr
\bold{4.} \eqn{P_i(f)} - itself: \eqn{PM_B * f + PN_B *(1-f) = AF*CN}, where \eqn{f=f_{snv}=f_{cnv}; PM_B<=pm}.\cr\cr

Under 1), SNV and CNV are completely independent as they are never co-propagated during the same clonal expansion. Under 2) and 3), SNV and CNV are partially dependent, yet present in two distinct subpopulations. Under 4), both the SNV and an CNV at \eqn{l} were propagated during the same clonal expansion.\cr 
}
\value{List with four components:
	  \item{p}{
	         The probability that the point mutation/CNV is present in a fraction \eqn{f} of cells, for each input frequency \eqn{f} in parameter \eqn{freq}.
    }         
	  \item{bestF}{
                 The cellular frequency that best explains the observed allele frequency and/or copy number.
    }
}
\references{Noemi Andor, Julie Harness, Sabine Mueller, Hans Werner Mewes and Claudia Petritsch. (2013) ExPANdS: Expanding Ploidy and Allele Frequency on Nested Subpopulations. Bioinformatics.}
\author{Noemi Andor}
\examples{
freq=seq(0.1,1.0,by=0.01);
cfd=cellfrequency_pdf(af=0.26,cnv=1.95,pnb=0,freq=freq, max_PM=6)
plot(freq,cfd$p,type="l",xlab="f",ylab="P(f)");
}
