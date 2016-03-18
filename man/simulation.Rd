\name{simulation}
\alias{simulation}
\docType{data}
\title{Simulated heterogeneous samples}
\description{A total of 50 samples with various numbers of subpopulations per sample were simulated at variable noise rates and constant number of 200 mutations per sample.}
\usage{data(simulation)}
\format{
        List 50 entries - one per simulated sample. Subpopulation composition can be predicted for each sample and the predictions compared to the the simulated entries:\cr
         \bold{snv} - the matrix of simulated point mutations (including ground truth columns SP*, PM*). \cr
         \bold{cbs} - the matrix of simulated point copy number segments (including ground truth columns SP*). \cr
         \bold{spstats} - matrix of subpopulation statistics (ground truth). 
}

\keyword{datasets}
