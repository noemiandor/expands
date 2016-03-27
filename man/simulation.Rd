\name{simulation}
\alias{simulation}
\docType{data}
\title{Simulated heterogeneous samples}
\description{A total of 50 samples with various numbers of subpopulations per sample were simulated at variable noise rates and constant number of 200 mutations per sample.}
\usage{data(simulation)}
\format{
        List with 50 entries - one per simulated sample. Subpopulation composition can be predicted for each sample and the predictions compared to the the simulated entries:\cr
         \bold{snv} - the matrix of simulated point mutations (including ground truth columns SP*, PM*). \cr
         \bold{cbs} - the matrix of simulated copy number segments (including ground truth columns SP*). \cr
         \bold{spstats} - matrix of subpopulation statistics (ground truth). 
}

\keyword{datasets}
\examples{
data(simulation)
snvcols=c("chr", "startpos", "CN_Estimate", "AF_Tumor","PN_B")
cbscols=c("chr", "startpos", "endpos")
sI=1:50;#set to 1:200 to run on entire simulation
#out=runExPANdS(simulation[[1]]$snv[sI,snvcols],simulation[[1]]$cbs[,cbscols],plotF = 0);
#truePhy=buildPhylo(simulation[[1]]$cbs,outF='truePhylo'); ##simulated 
#predPhy=buildPhylo(out$ploidy,outF='truePhylo'); ##predicted
#par(mfrow=c(1,2))
#plot(truePhy$tree,cex=2,main='simulated')
#plot(predPhy$tree,cex=2,main='predicted')
}