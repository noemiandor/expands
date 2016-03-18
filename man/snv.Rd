\name{snv}
\alias{snv}
\docType{data}
\title{Single Nucleotide Variations}
\description{Somatic mutations and Loss of Heterozygosity (LOH) of a Glioblastoma tumor (TCGA-06-0152-01).}
\usage{data(snv)}
\format{
	Numeric matrix with 773 rows (one per mutation) and 7 columns:\cr
	 \bold{chr} - the chromosome \cr
	 \bold{startpos} - genomic position \cr
	 \bold{endpos} - same as above \cr
	 \bold{REF} - ASCII code of the reference nucleotide (in hg18/hg19) \cr
	 \bold{ALT} - ASCII code of the B-allele nucleotide \cr
	 \bold{AF_Tumor} - allele frequeny of B-allele \cr
	 \bold{PN_B} - ploidy of B-allele in normal cells. A value of 0 indicates that the mutation has only been detected in the tumor sample (i.e. somatic mutations). A value of 1 indicates that the mutation is also present in the normal (control) sample, albeit at reduced allele frequency (i.e. mutation is consequence of LOH). Other mutations should not be included.
}

\source{Data derived from The Cancer Genome Atlas (TCGA).}


\keyword{datasets}
