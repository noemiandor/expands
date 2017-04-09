\name{buildPhylo}
\alias{buildPhylo}
\title{Relations between subpopulations}
\description{Predicts phylogenetic relations between subpopulations from subpopulation specific copy number profiles.}
\usage{buildPhylo(ploidy, outF, treeAlgorithm="bionjs", dm=NA, add="Germline")
}\arguments{
  \item{ploidy}{
       Ploidy-matrix in which each row corrsponds to a copy number segment. Has to contain at least one column for each predicted subpopulation. Subpopulation columnnames must be labeled SP_xx, where xx is the size of the corresponding subpopulation. Ploidy-matrix can be obtained by calling \code{\link{assignQuantityToSP}}. 
}
  \item{outF}{
       Prefix of file to which phylogeny will be saved.
}
  \item{treeAlgorithm}{
       Neighbor joining algorithm used for phylogeny reconstruction (from library ape). Options: bionjs (default), njs.
}
  \item{dm}{
       Optional matrix in which each row corresponds to a mutation. Only mutations located on autosomes should be included. Columns in dm must be labeled and must include:\cr
        \bold{SP} - the subpopulation to which the point mutation has been asssigned.\cr
        \bold{SP_cnv} - the subpopulation to which the CNV (overlapping with the point mutation) has been asssigned (if an CNV is present).\cr
        \bold{chr} - the chromosome on which each point mutation is located;\cr
        \bold{startpos} - the genomic position of each point mutation; \cr	
        \bold{PM} - the total ploidy of all alleles at the mutated genomic locus, in the assigned subpopulation.\cr
		    \bold{PM_B} - the ploidy of the B-allele at the mutated genomic locus, in the assigned subpopulation.\cr  
         If dm is available, an attempt will be made to assign every mutation to >1 subpopulation according to the inferred phylogenetic relations between subpopulations.
}
  \item{add}{
       Artificial subpopulation to be included in phylogeny (options: 'Germline','Consensus', NULL).
}
}
\details{
Reconstructs phylogenetic relationships between subpopulations using neighbor-joining algorithms provided by R-package 'ape'. Pairwise distances between subpopulations are calculated as the number of copy number segments for which both subpopulations have the same copy number, divided by the total number of copy number segments for which both subpopulations have available copy number information. Subpopulations with insufficient ploidy information are excluded from phylogeny.
}
\value{List with two fields:
        \item{tree}{An object of class "phylo" (library ape).}
        \item{dm}{The input matrix with each row representing a point mutation and additional columns:\cr
        	\bold{SP_xx} - where xx is the size of the corresponding subpopulation. Column entries contain a binary indicator of whether or not the point mutation in this row is present in SP_xx.
        	}
}

\author{Noemi Andor}


\seealso{
\code{\link{assignQuantityToSP}}
}

