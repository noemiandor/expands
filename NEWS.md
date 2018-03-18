		CHANGES IN EXPANDS VERSION 2.1.1

BUG FIXES

    o Fixed bug in buildPhylo that prevented recognition of subpopulations within the phylogenetic tree, thereby limiting assignment of any given mutation to no more than one subpopulation.



		CHANGES IN EXPANDS VERSION 2.1.0

BUG FIXES

    o Fixed undesired behaviour in Subpopulations.java where integer division was rounded thereby rejecting valid mutated and wild type copy number configurations.
    o Fixed bug in runExPANdS.R where the number of mutations per subpopulation was recorded incorrectly in .spstats output.
    o Fixed axis-range bug in plotSPs.R.

NEW FEATURES

    o Included parameter ploidy in runExPANdS.R, allowing specification of non-diploid background states (e.g. for near-triploid cell lines).
    o Parallel computing option for runExPANdS.R (parameter nc).

OTHER CHANGES	
    o Copy-neutral LOH are no longer modelled.
    o Fit of copy number status has a higher weight (x^2) than fit of point mutation status (x^1).
    o Likelihood of mutation events is no longer higher in cancer-cell populations known to be non-founders.
    o Binary option 'ambig' in assignQuantityToSP.R replaced with parameter 'e', allowing user to finetune amount of accepted ambiguity within a  continuous range (e==0 is equivalent to former ambig==F, while e>0 corresponds to ambig==T).



		CHANGES IN EXPANDS VERSION 2.0.0

BUG FIXES

    o Fixed bug in plotSPs(..): adjusted y-axis to include 0 (i.e. deletion scenario).

NEW FEATURES

    o Copy-neutral LOH are now modelled. 

OTHER CHANGES	

    o More robust calculation of cell frequency probabilities from kernel density estimates instead of Gaussian mixtures.
    o Complete remodelling of java core to better transparency of the four different evolutionary scenarios.


	
		CHANGES IN EXPANDS VERSION 1.7.2		
BUG FIXES

    o Fixed bug in function plotSPs: range of %maxP values was too large to be covered by the gray color map, resulting in subpopulations not being displayed. %maxP values are now scaled between 0 and 50.



                CHANGES IN EXPANDS VERSION 1.7.1
BUG FIXES

    o Summary file of detected subpopulations was not included to the list of output files.



		CHANGES IN EXPANDS VERSION 1.7


NEW FEATURES

    o New subpopulation plotting options (function plotSPs).
    o New phylogeny plotting option: user can choose between consensus and germline population to be included as control.
    o Two additional output files created: subpopulation specific ploidy matrix and a summary file of all detected subpopulations, including information on the ancestor and the closest descendant of each subpopulation.
    o 50 simulated samples of various genetic complexities and noise levels included for user testing.
    o Whenever possible, SNVs are assigned not only to the subpopulation in which they first occurred, but also to descending subpopulations.

BUG FIXES

    o Fixed bug where rare scenarios would lead to wrong subpopulation specific ploidy assignments at specific mutated loci. Once the subpopulation size was assigned to a mutation, EXPANDS was using the magnitude of minuscule fluctuations in cellular frequencies to choose the “winner” among alternative ploidy scenarios (associated with the chosen subpopulation size). Instead it should be using the magnitude of the residual from equation 3-solutions to choose the ploidy scenario that best explains the data.
 

OTHER CHANGES

    o Improved subpopulation detection specificity, in particular for small subpopulations, by removing “doublet” subpopulations - where the heterozygous and homozygous state of a single subpopulation was wrongly split into two subpopulations.
    o Improved mutation assignment to subpopulations. Within individual mutations, the sum of probabilities across input-cellfrequencies was always normalized to 1. This is important for clustering and used to be fine for the mutation assignment step as well for as long as there was just one distribution per mutation. But ever since version 1.6, EXPANDS covers the different evolutionary scenarios (SNV simultaneous with/before/after CNV, etc.) – so there are up to 4 distributions per mutation. These are no longer comparable amongst each other when normalizing. Normalization to 1 needs to be skipped for the mutation assignment step.



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
