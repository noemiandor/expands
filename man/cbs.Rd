\name{cbs}
\alias{cbs}
\docType{data}
\title{Matrix of copy number fragments}
\description{Copy number segments as obtained by circular binary segmentation. Data is derived from a Glioblastoma tumor (TCGA-06-0152-01).}
\usage{data(cbs)}
\format{
	Numeric matrix with 120 rows (one per copy-number segment) and 4 columns:\cr
	 \bold{chr} - the chromosome \cr
	 \bold{startpos} - genomic position at which copy-number segment starts.\cr
	 \bold{endpos} - genomic position at which copy-number segment ends.\cr
	 \bold{CN_Estimate} - average copy-number of the segment among all cells.
}

\source{Data derived from The Cancer Genome Atlas (TCGA).}


\keyword{datasets}
