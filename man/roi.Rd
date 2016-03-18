\name{roi}
\alias{roi}
\docType{data}
\title{Regions of interest}
\description{For internal use only. Default regional boundary for mutations included during clustering, comprising ca. 468 MB centered on the human exome. Relevant if number of input mutations exceeds user defined threshold (often applies to whole genome sequencing data). A saved image of this object is in sysdata.rda.}
\format{
	Numeric matrix in which each row corresponds to a genomic segment. \cr
	Columns: \cr
        \bold{chr} - the chromosome of the segment ;\cr
        \bold{start} - the first genomic position of the segment;\cr
        \bold{end} - the last genomic position of the segment.

}

\source{Data derived from human SureSelectExome_hg19 50MB library kit annotation.}


\seealso{
\code{\link{runExPANdS}}
}

\keyword{datasets}
