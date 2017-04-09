\name{assignQuantityToSP}
\alias{assignQuantityToSP}
\title{Quantity assignment (ploidy) to subpopulations}
\description{Assigns quantities to predicted subpopulations. Currently, the only assignable quantity are subpopulation specific ploidies for the input genome segments.
}
\usage{assignQuantityToSP(cbs, dm, C=list(sps=c("SP","SP_cnv"),pms = c("PM","PM_cnv")), ambig=F)}
\arguments{
  \item{cbs}{
	Matrix in which each row corrsponds to a copy number fragment as computed by a circular binary segmentation algorithm. Has to contain at least the following columnnames:\cr
	\bold{chr} - chromosome;\cr
	\bold{startpos} - the first genomic position of a copy number segment; \cr
	\bold{endpos} - the last genomic position of a copy number segment; \cr
	\bold{CN_Estimate} - the copy number estimated for each segment (weighted average value across all subpopulations in the sample). 
}
  \item{dm}{
	Matrix in which each row corrsponds to a mutation. Has to contain at least the following columnnames: \cr 
	\bold{chr} - the chromosome on which each point mutation is located;\cr
	\bold{startpos} - the genomic position of each mutation;\cr
	\bold{SP} - the subpopulation to which the point mutation has been assigned;\cr
	\bold{SP_cnv} - the subpopulation with a copy number variation within the same genomic segment in which SP has a point mutation;\cr
	\bold{PM} - ploidy of the subpopulation with the point mutation (SP);\cr
	\bold{PM_cnv} - ploidy of the subpopulation with the copy number variation (SP_cnv).
}
  \item{C}{List referencing column names in the mutation matrix, with two fields:\cr
  	\bold{sps} - the column names holding subpopulation sizes (typically "SP", "SP_cnv"); \cr
  	\bold{pms} - the column names holding the ploidies assigned for each subpopulation (typically "PM", "PM_cnv"). \cr
}
\item{ambig}{
        Whether or not to assign ploidy to a subpopulation, \eqn{SP_i}, for a segment containing multiple \eqn{SP_i} specific ploidies, at least two of which have distinct magnitudes. If set to TRUE, the median ploidy is assigned as segment ploidy. Setting this parameter to TRUE is not recommended as the output will include segment-assignments where subpopulation specific ploidy is ambiguous. Repeating circular binary segmentation with less stringent parameters typically reduces segment length and thus the prevalence of ambiguous assignments.
}
}

\value{The input copy number matrix with one additional column for each predicted subpopulation:\cr
	\bold{SP_xx} - where xx is the size of the corresponding subpopulation. Column entries contain the ploidy of each segment in SP; Value <NA> indicates that no ploidy could be inferred for the segment in this subpopulation (either because the subpopulation had no point mutations/CNVs within the segment, or because it had multiple, ambiguous ploidy assignments within the segment).
}

\author{Noemi Andor}
