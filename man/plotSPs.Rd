\name{plotSPs}
\alias{plotSPs}
\title{Subpopulation Visualization}
\description{Plots coexistent subpopulations determined by ExPANdS.}
\usage{
plotSPs(dm, sampleID=NA, cex=0.5, legend="CN_Estimate", orderBy="chr", rawAF=F)
}
\arguments{
  \item{dm}{
	 Matrix in which each row corresponds to a point mutation (for example, the matrix output of \code{\link{assignMutations}}). Has to contain at least the following column names:\cr
        \bold{chr} - the chromosome on which each mutation is located; \cr
        \bold{startpos} - the genomic position of each mutation;\cr
        \bold{AF_Tumor} - the allele-frequency of each mutation;\cr
        \bold{CN_Estimate} - the absolute copy number estimated for each segment;\cr
        \bold{PN_B} - the count of the B-allele in normal cells (binary variable: 1 if the mutation is a germline variant, 0 if somatic);\cr
	      \bold{SP} - the subpopulation to which each point mutation has been assigned (as fraction of cells in the tumor bulk);\cr
	      \bold{\%maxP} -	the confidence with which the mutation has been assigned to the corresponding subpopulation;\cr
	      \bold{SP_cnv} - the subpopulation to which the CNV has been assigned; \cr
	      \bold{PM} - the total count of all alleles at the mutated genomic locus, in subpopulation SP.\cr
	      \bold{PM_cnv} - the total count of all alleles at the mutated genomic locus in subpopulation SP_cnv.\cr
	      \bold{PM_B} - the count of the mutated allele in subpopulation SP.\cr
	      \bold{PM_B_cnv} - the count of the mutated allele in subpopulation SP_cnv.
}
  \item{sampleID}{
	The name of the sample in which the mutations have been detected.
}
  \item{cex}{
	The amount by which plotting text and symbols should be magnified relative to the default. See also \code{help(par).}
}

  \item{legend}{
	Allele frequencies and subpopulation specific copy numbers are colored based on the chromosome on which the mutation is located (option: 'chr') or based on the average copy number of the locus in the sample (option: 'CN_Estimate').
}        
  \item{orderBy}{
	Loci within a subpopulation are sorted by genomic location (option 'chr') or by the confidence with which they have been assigned to the subpopulation (option '\%maxP').
}
  \item{rawAF}{
	Specifies whether the allele frequency of SNVs should be adjusted relative to the assigned subpopulation (options: true, false).
}
}

\value{
For each point mutation (x-axis) the function displays:\cr
     - the size of the subpopulation to which the mutation has been assigned (squares). Each square is colored based on the confidence with which the mutation has been assigned to the corresponding subpopulation (black - highest, white - lowest).\cr
     - the total count of all alleles at the mutated genomic locus in that subpopulation (dots).\cr
     - only for loci with an CNV and an SNV each in distinct subpopulations: the total count of all alleles at the mutated genomic locus in the subpopulation which harbors the CNV (crosses).\cr
     - the allele frequency of the mutation (stars - somatic mutations, triangles - loss of heterozygosity).
}

\author{Noemi Andor}
