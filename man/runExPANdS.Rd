\name{runExPANdS}
\alias{runExPANdS}
\title{Main Function}
\description{Given a set of mutations, ExPANdS predicts the number of clonal expansions in a tumor, the size of the resulting subpopulations in the tumor bulk and which mutations accumulate in a cell prior to its clonal expansion. Input-parameters SNV and CBS hold the paths to tab-delimited files containing the point mutations and the copy numbers respectively. Alternatively SNV and CBS can be read into the workspace and passed to \code{runExPANdS}  as numeric matrices. The robustness of the subpopulation predictions by ExPANdS increases with the number of mutations provided. It is recommended that SNV contains at least 200 point mutations to obtain stable results.}
\usage{runExPANdS(SNV, CBS, maxS=0.7, max_PM=6, min_CF=0.1, p=NA, ploidy=2, 
                  nc=1, plotF=2, snvF=NULL, maxN=8000, region=NA, verbose=T)}
\arguments{
  \item{SNV}{
	Matrix in which each row corresponds to a point mutation. Only mutations located on autosomes should be included. Columns in SNV must be labeled and must include: \cr
	\bold{chr} - the chromosome on which each mutation is located;\cr
	\bold{startpos} - the genomic position of each mutation; \cr
	\bold{AF_Tumor} - the allele-frequency of each mutation;\cr
	\bold{PN_B} - count of B-allele in normal cells. A value of 0 indicates that the variant has only been detected in the tumor sample (i.e. somatic mutation). A value of 1 indicates that the variant is also present in the normal (control) sample, albeit at reduced allele frequency (i.e. this is a germline variant, which passed the calling filter due to the presence of an LOH event). Mutations, for which the allele frequency in the tumor sample is lower than the corresponding allele frequency in the normal sample, should not be included.
}
  \item{CBS}{
	Matrix in which each row corresponds to a copy number segment. CBS is typically the output of a circular binary segmentation algorithm. Columns in CBS must be labeled and must include:\cr
	\bold{chr} - chromosome; \cr
	\bold{startpos} - the first genomic position of a copy number segment;\cr 
	\bold{endpos} - the last genomic position of a copy number segment; \cr
	\bold{CN_Estimate} - the absolute copy number estimated for each segment.
}
  \item{maxS}{
	Upper threshold for the noise score of subpopulation detection. Only subpopulations identified at a score below \eqn{maxS} are kept.
}
  \item{max_PM}{
	Upper threshold for the number of amplicons per mutated cell. Increasing the value of this variable is not recommended unless extensive depth and breadth of coverage underlie the measurements of copy numbers and allele frequencies. See also \code{\link{cellfrequency_pdf}.}
}
  \item{min_CF}{
	Lower boundary for the cellular prevalence interval of a mutated cell. Mutations for which allele frequency * copy number are below \eqn{min_{CellFreq}}, are excluded from further computation. Decreasing the value of this variable is not recommended unless extensive depth and breadth of coverage underlie the measurements of copy numbers and allele frequencies.}
  \item{p}{
	Precision with which subpopulation size is predicted, a small value reflects a high resolution and can lead to a higher number of predicted subpopulations.}
  \item{plotF}{
	Option for displaying a visual representation of the identified subpopulations (0 - no display; 1 - display subpopulation size; 2 - display subpopulation size and phylogeny).}
  \item{snvF}{
	Prefix of file to which predicted subpopulation composition will be saved. Default: the name of the file from which mutations have been read or "out.expands" if input mutations are not handed over as file path.}
  \item{maxN}{
	Upper limit for number of point mutations used during clustering. If number of user supplied point mutations exceeds \eqn{maxN}, the clustering of cellular frequency distributions will be restricted to point mutations found within \eqn{region}.}
  \item{region}{
	Regional boundary for mutations included during clustering. 
	Matrix in which each row corresponds to a genomic segment. Columns must include: \cr
	  \bold{chr} - the chromosome of the segment;\cr
	  \bold{start} - the first genomic position of the segment;\cr
	  \bold{end} - the last genomic position of the segment.\cr
	Default: SureSelectExome_hg19, comprising ca. 468 MB centered on the human exome. Alternative user supplied regions should also be coding regions, as the selective pressure is higher as compared to non-coding regions.
}
  \item{ploidy}{
        The background ploidy of the sequenced sample (default: 2). Changing the value of this parameter is not recommended. Dealing with cell lines or tumor biopsies of very   high (>=0.95) tumor purity is a necessary but not sufficient condition to change the value of this parameter.}
  \item{nc}{
        The number of nodes to be forked to run R in parallel.} 
  \item{verbose}{
        Give a more verbose output.}
}
\value{List with fields:
        \item{finalSPs}{Matrix of predicted subpopulations. Each row corresponds to a subpopulation and each column contains information about that subpopulation, such as the size in the sequenced tumor bulk (column \bold{Mean Weighted}) and the noise score at which the subpopulation has been detected (column \bold{score}).}
        \item{dm}{Matrix containing the input mutations with at least seven additional columns:\cr
	\bold{SP} - the subpopulation to which the point mutation has been assigned; \cr
	\bold{SP_cnv} - the subpopulation to which the CNV has been assigned (if an CNV exists at this locus); \cr
	\bold{\%maxP} - the confidence of mutation assignment.\cr
	\bold{f} - Deprecated. The maximum likelihood cellular prevalence of this point mutation, before it has been assigned to SP. This value is based on the copy number and allele frequency of the mutation exclusively and is independent of other point mutations. Column SP is less sensitive to noise and considered the more accurate estimation of cellular mutation prevalence.\cr
	\bold{PM} - the total count of all alleles in the subpopulation harboring the point mutation (SP).\cr
	\bold{PM_B} - the count of the B-allele in the subpopulation harboring the point mutation (SP).\cr	     
	\bold{PM_cnv} - the total count of all alleles in the subpopulation harboring an CNV (SP_cnv).\cr
	\bold{PM_B_cnv} - the count of the B-allele, in the CNV harboring subpopulation (SP_cnv).\cr
	If phylogeny reconstruction was successful, matrix includes one additional column for each subpopulation from the phylogeny, indicating whether or not the point mutation is present in the corresponding subpopulation.\cr}
	\item{densities}{Matrix as obtained by \code{\link{computeCellFrequencyDistributions}.} Each row corresponds to a mutation and each column corresponds to a cellular frequency. Each value \eqn{densities[i,j]} represents the probability that mutation \eqn{i} is present in a fraction \eqn{f} of cells, where \eqn{f} is given by: \eqn{colnames(densities[,j]).}}
	\item{sp_cbs}{Matrix as obtained by \code{\link{assignQuantityToSP}.} Each row corresponds to a copy number segment, e.g. as obtained from a circular binary segmentation algorithm. Includes one additional column for each predicted subpopulation, containing the copy number of each segment in the corresponding subpopulation.}
	\item{tree}{An object of class "phylo" (library ape) as obtained by \code{\link{buildPhylo}.} Contains the inferred phylogenetic relationships between subpopulations.}
}
\references{Noemi Andor, Julie Harness, Sabine Mueller, Hans Werner Mewes and Claudia Petritsch. (2013) ExPANdS: Expanding Ploidy and Allele Frequency on Nested Subpopulations. Bioinformatics.}
\author{Noemi Andor}
\examples{
data(snv);
data(cbs);
maxS=2.5;
set.seed(4); idx=sample(1:nrow(snv), 60, replace=FALSE);
#out= runExPANdS(snv[idx,], cbs, maxS);
}
