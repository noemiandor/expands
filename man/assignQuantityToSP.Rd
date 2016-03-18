\name{assignQuantityToSP}
\alias{assignQuantityToSP}
\title{Quantity assignment (ploidy) to subpopulations}
\description{Assigns quantities to predicted subpopulations. Currently, the only assignable quantity are subpopulation specific ploidies for the input genome segments (obtained for example by circular binary segmentation).
}
\usage{assignQuantityToSP(cbs, dm, colName = "PM_cnv", keepAmbigSeg=FALSE)}
\arguments{
  \item{cbs}{
	Matrix in which each row corrsponds to a copy number fragment as computed by a circular binary segmentation algorithm. Has to contain at least the following columnnames:\cr
	\bold{chr} - chromosome;\cr
	\bold{startpos} - the first genomic position of a copy number segment; \cr
	\bold{endpos} - the last genomic position of a copy number segment; \cr
	\bold{CN_Estimate} - the copy number estimated for each segment (average value across all subpopulations in the sample). 
}
  \item{dm}{
	Matrix in which each row corrsponds to a mutation. Has to contain at least the following columnnames: \cr 
	\bold{chr} - the chromosome on which each point mutation is located;\cr
	\bold{startpos} - the genomic position of each mutation;\cr
	\bold{SP} - the subpopulation to which the point mutation has been assigned;\cr
	\bold{SP_cnv} - the subpopulation with copy number variation (CNV) in the same genomic segment in which SP has a point mutation;\cr
	\bold{PM} - ploidy of the subpopulation with the point mutation (SP);\cr
	\bold{PM_cnv} - ploidy of genomic segment (overlapping the point mutation), in the subpopulation with CNV (SP_cnv).
}
  \item{colName}{
	The subpopulation specific value assigned to each copy number fragment. Possible values: PM_cnv, PM, PM_B. Default: PM_cnv.
}
\item{keepAmbigSeg}{
        Whether or not to assign ploidy to a subpopulation, \eqn{SP_i}, for a segment containing multiple \eqn{SP_i} specific ploidities, at least two of which have distinct magnitudes. If set to TRUE, the median ploidy is assigned as segment ploidy. Setting this parameter to TRUE is not recommended as the output will include segment-assignements where subpopulation specific ploidy is ambiguous.Recommend repeating circular binary segmentation with less stringent parameters instead, to reduce segment length and thus the prevalence of ambiguous assignements. Default: FALSE.
}
}

\value{List with two fields:
        \item{dm}{The input mutation matrix with missing values in column PM imputed from information in copy number matrix when possible.}
	\item{ploidy}{The input copy number matrix with one additional column for each predicted subpopulation:\cr
	\bold{SP_xx} - where xx is the size of the corresponding subpopulation (SP). Column entries contain the ploidy of each segment in SP; Value <NA> indicates that no ploidy could be inferred for the segment in this subpopulation (either because SP had no point mutations/CNVs within the segment, or because SP had multiple, ambiguous ploidy assignments within the segment).
}
}

\author{Noemi Andor}
