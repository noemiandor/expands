\name{buildMultiSamplePhylo}
\alias{buildMultiSamplePhylo}
\title{Relations between inter- and intra-sample subpopulations}
\description{Predicts phylogenetic relations between subpopulations from subpopulation specific copy number and point mutation profiles, while including information about sample origin of each subpopulation. This function differs from \code{\link{buildPhylo}} in that it integrates the subpoulations predicted in multiple, geographically distinct tumor-samples into one common phylogeny and in that it includes point mutations in addition to copy number variations to infer inter-sample phylogenetic relations.}
\usage{buildMultiSamplePhylo(samGr,out,treeAlgorithm="bionjs",ambigSg=F,plotF=1,spRes=1)}\arguments{
  \item{samGr}{
    List with three fields:\cr
	\bold{cbs} - Matrix in which each row corresponds to a copy number segment. CBS is typically the output of a circular binary segmentation algorithm. Columns in CBS must be labeled and must include chr, startpos, endpos and CN_Estimate (see \code{\link{cbs}}).\cr
	\bold{sps} - Output of \code{\link{runExPANdS}}. Matrix in which each row corresponds to a somatic mutations. Columns must include: chr - the chromosome on which each mutation is located; startpos - the genomic position of each mutation; SP - the subpopulation to which the mutation has been asssigned; PM - the total ploidy of all alleles at the mutated genomic locus, in the assigned subpopulation; PM_B - the ploidy of the B-allele at the mutated genomic locus, in the assigned subpopulation; CN_Estimate - the average ploidy (among all cells) of the locus in which the mutation is embedded (see also \code{\link{assignQuantityToMutation}}).\cr  
	\bold{labels} - Label denoting sample origin of each subpopulation matrix. 
    Entry is mandatory for each geographical sample.\cr
  }
  \item{out}{
    Prefix of file to which multi-sample phylogeny will be saved.
  }
  \item{treeAlgorithm}{
    Neighbor joining algorithm used for phylogeny reconstruction (from library ape). Options: bionjs (default), njs.
  }
  \item{ambigSg}{
    Input parameter "keepAmbigSeg" for called function: \code{\link{assignQuantityToSP}}. 
  }
  \item{plotF}{
    Option for displaying a visual representation of the phylogenetic tree (0 - no display; 1 - display). Default: 1.
  }
  \item{spRes}{
    Option on whether or not to ignore the subpopulations calculated for each sample and instead treat every geographical tumor-sample as one single tumor-metapopulation (Default value: 1 - subpopulation resolution; 0 - metapopulation resolution). 
  }
}
\details{
	Reconstructs phylogenetic relationships between subpopulations using neighbor-joining algorithms provided by R-package 'ape'. Pairwise distances between subpopulations \eqn{i} and \eqn{j} are calculated as: \cr
	\eqn{d_{ij}:=(cnv_{i=j} + snv_{i=j})/(cnv_{ij}+snv_{ij})}, where \eqn{ cnv_{i=j} } is the number of copy number segments for which subpopulations \eqn{i} and \eqn{j} have the same copy number; \eqn{ snv_{i=j} } is the number of point mutations for which subpopulations \eqn{i} and \eqn{j} have the same mutation status and  \eqn{ cnv_{ij}, snv_{ij} } are the total number of copy number segments and mutations respectively, for which both subpopulations have available information. Subpopulations with insufficient ploidy and point mutations information are excluded from phylogeny. 
}
\value{An object of class "phylo" (library ape).
}

\author{Noemi Andor}


\seealso{
  \code{\link{buildPhylo}}
}
