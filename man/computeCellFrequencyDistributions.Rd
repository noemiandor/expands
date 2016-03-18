\name{computeCellFrequencyDistributions}
\alias{computeCellFrequencyDistributions}
\title{Gathering of cell frequency probability distributions}
\description{Computes the probability distributions of cell frequencies, by calling \code{\link{cellfrequency_pdf}} for each mutation separately.}
\usage{computeCellFrequencyDistributions(dm, max_PM=6, precision, min_CellFreq=0.1)}
\arguments{
  \item{dm}{
	 Matrix in which each row corrsponds to a mutation. Has to contain at least the following columnnames: \cr 
	\bold{chr} - the chromosome on which each mutation is located; \cr
	\bold{startpos} - the position of each mutation;\cr 
	\bold{AF_Tumor} - the allele-frequency of each mutation;\cr
	\bold{PN_B} - the ploidy of the B-allele in normal cells (binary variable: 1 if the mutation is a germline variant, 0 if somatic).
}
  \item{max_PM}{
        Upper threshold for the number of amplicons per mutated cell (default: 6). See also \code{\link{cellfrequency_pdf}.} 
}
  \item{precision}{
	Precision with which subpopulation size is predicted, a small value reflects a high resolution and can lead to a higher number of predicted subpopulations.
}
  \item{min_CellFreq}{
        Lower boundary for the prevalence of a mutated cell (default: 0.1).}
}

\value{List with three fields:
	\item{freq}{The cellular frequencies for which probabilities are computed.}
	\item{densities}{Matrix in which each row corresponds to a point mutation and each column corresponds to a cellular frequency. Each value \eqn{densities[i,j]} represents the probability that mutation \eqn{i} is present in a fraction \eqn{freq[j]} of cells.
}
	\item{dm}{The input matrix with column \eqn{f} updated according to the cellular frequency that best explains the observed allele frequency and copy number.}
}

\author{Noemi Andor}


