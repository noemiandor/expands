\name{assignQuantityToMutation}
\alias{assignQuantityToMutation}
\title{Quantity assignment (copy number) to mutations}
\description{Assigns a quantity to each mutated locus. Currently, the only assignable quantity is the average copy number (among all cells) of the locus in which the mutation is embedded. 
}
\usage{assignQuantityToMutation(dm, cbs, quantityColumnLabel="CN_Estimate", verbose = T)}
\arguments{
  \item{dm}{
	Matrix in which each row corresponds to a mutation. Has to contain at least the following column names: \cr 
	\bold{chr} - the chromosome on which each mutation is located;\cr
	\bold{startpos} - the genomic position of each mutation.
}
  \item{cbs}{
	Matrix in which each row corresponds to a copy number segment as calculated by a circular binary segmentation algorithm. Has to contain at least the following column names:\cr
	\bold{chr} - chromosome;\cr
	\bold{startpos} - the first genomic position of a copy number segment; \cr
	\bold{endpos} - the last genomic position of a copy number segment; \cr
	\bold{CN_Estimate} - the copy number estimated for each segment. 
}
  \item{quantityColumnLabel}{
	The name of the new column. Valid options are: \code{FPKM, CN_Estimate.} 
}
\item{verbose}{
        Give a more verbose output.
}
}

\value{
\item{dm}{
       The input matrix with three additional columns:\cr
	\bold{quantityID} - the ID of the assigned quantity; \cr
	\bold{quantityColumnLabel} - the quantity.
}
}

\author{Noemi Andor}
\examples{
data(cbs)
data(snv)
dm=assignQuantityToMutation(snv,cbs,quantityColumnLabel="CN_Estimate")
}
