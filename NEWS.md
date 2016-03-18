		CHANGES IN EXPANDS VERSION 1.6.1

BUG FIXES

    o Fixed bug in assignQuantityToSP.R to take into account that column "PM" in .sps matrix can be NA
    o Fixed bug in buildMultiSamplePhylo.R to ensure merged .sps  matices always have consistent column dimensions

OTHER CHANGES

    o Function runExPANdS no longer substitutes existing dots (“.”) with underscores (“_”) for the .tree and .dist output files 



		CHANGES IN EXPANDS VERSION 1.6


NEW FEATURES

    o So far mutations had been assigned to maximal one subpopulation. However mutations may not be exclusive to the assigned subpopulation but may also be present in smaller, descending subpopulations. This version of expands decides whether or not this is the case leveraging the predicted phylogenetic structure of the subpopulation composition. 
    o Included homozygous deletion as potential scenario when modeling (SNV,CNV) pairs that have overlapping genomic location, but are propagated during distinct clonal expansions. 

BUG FIXES

    o Incorrect data matrix conversion could occur when handing over a non-numerical file as parameter "SNV" to function runExPANdS. This is now fixed.
    o Fixed bug of occasionally having cellular frequencies >1 as output of subpopulation sizes. This bug fix was the consequence of optimizing the solution through which sensitivity at cell-frequency distribution margins is improved. Need for improvement was because subpopulation detection sensitivity correlates to centrality of subpopulation size during clustering. Tolerance of copy number and allele frequency measurement errors must be higher for marginal cell-frequencies than for central cell-frequencies, in order to counteract the reduced cluster detection sensitivity at the cell-frequency distribution margins. This is only relevant during subpopulation detection (clustering), uniform error tolerance still applies during SNV assignment.



		CHANGES IN EXPANDS VERSION 1.5


NEW FEATURES

    o Additional optional parameter "min_CellFreq" provided for function runExPANdS, specifying the minimum cellular prevalence of mutations to be included for subpopulation predictions.
    o Additional function "buildMultiSamplePhylo" available, which integrates the subpopulations predicted in multiple, geographically distinct tumor samples into one common phylogeny.


OTHER CHANGES

    o Filtered loci with high-level amplifications, according to max_PM setting. This reduces unnecessary processing time, as assignment of mutations within amplified regions to subpopulations is not successful.
