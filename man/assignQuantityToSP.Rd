\name{assignQuantityToSP}
\alias{assignQuantityToSP}
\title{Quantity assignment (copy number) to subpopulations}
\description{Assigns quantities to predicted subpopulations. Currently, the only assignable quantity are subpopulation specific copy number states for the input genome segments.
}
\usage{assignQuantityToSP(cbs, dm, C=list(sp=c("SP","SP_cnv"), pm=c("PM","PM_cnv")), e=1, v=T)}
\arguments{
  \item{cbs}{
	Matrix in which each row corresponds to a copy number fragment as computed by a circular binary segmentation algorithm. Has to contain at least the following column names:\cr
	\bold{chr} - chromosome;\cr
	\bold{startpos} - the first genomic position of a copy number segment; \cr
	\bold{endpos} - the last genomic position of a copy number segment; \cr
	\bold{CN_Estimate} - the copy number estimated for each segment (weighted average value across all subpopulations in the sample). 
}
  \item{dm}{
	Matrix in which each row corresponds to a mutation. Has to contain at least the following column names: \cr 
	\bold{chr} - chromosome on which each point mutation is located;\cr
	\bold{startpos} - genomic position of each mutation;\cr
	\bold{SP} - subpopulation to which the point mutation has been assigned;\cr
	\bold{SP_cnv} - subpopulation with a copy number variation within the same genomic segment in which SP has a point mutation;\cr
	\bold{PM} - total count of all alleles in the subpopulation with the point mutation (SP);\cr
	\bold{PM_cnv} - total count of all alleles in the subpopulation with the copy number variation (SP_cnv).
}
  \item{C}{List referencing column names in the mutation matrix, with two fields:\cr
  	\bold{sp} - column names holding subpopulation sizes (typically "SP", "SP_cnv"); \cr
  	\bold{pm} - column names holding the total allele counts assigned for each subpopulation (typically "PM", "PM_cnv"). \cr
}
\item{e}{
        Maximum variance of subpopulation specific copy numbers for a given segment, above which segment will remain unassigned for the corresponding subpopulation. Determines whether or not to assign copy number to a subpopulation, \eqn{SP_i}, for a segment containing multiple \eqn{SP_i} specific copy numbers, at least two of which are distinct.
}
\item{v}{
        Give a more verbose output.
}
}

\value{The input copy number matrix with one additional column for each predicted subpopulation:\cr
	\bold{SP_xx} - where xx is the size of the corresponding subpopulation. Column entries contain the copy number of each segment in SP; Value <NA> indicates that no copy number could be inferred for the segment in this subpopulation (either because the subpopulation had no point mutations/CNVs within the segment, or because it had multiple, ambiguous copy number assignments within the segment).
}

\author{Noemi Andor}
