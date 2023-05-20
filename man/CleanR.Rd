\name{CleanR}
\alias{CleanR}
\title{Cluster enhancement for testing intermodal associations}
\usage{
CleanR(
  ymat,
  xmat,
  distmat,
  sacf = "exp",
  max.radius = 20,
  nperm = 5000,
  alpha = 0.05,
  alternative = c("two.sided", "less", "greater"),
  seed = NULL,
  partition = T,
  npartition = NULL,
  parallel = F,
  ncores = 1
)
}

\arguments{
\item{ymat}{A V times N matrix (V: # of vertices, N: # of images) for the first modality.}

\item{xmat}{A V times N matrix (V: # of vertices, N: # of images) for the second modality. \code{xmat} and \code{ymat} should have their subject IDs aligned by row.}


\item{mod}{A N times p covariate matrix (p: # of covariates). It can be generated easily by using the \code{model.matrix()} function.}

\item{distmat}{A V times V distance matrix}

\item{sacf}{A spatial autocorrelation function. The exponential function (\code{exp}) is assumed as a default. Other choices include \code{gau} (Gaussian) or \code{mix} (mixture of exponential and Gaussian).}

\item{max.radius}{The maximum radius for cluster enhancement. 20 is assumed as a default.}

\item{nperm}{The number of permutations to be used. At least 5000 permutation is recommended. Default: 5000}

\item{alpha}{A desired FWER. Default: 0.05.}

\item{alternative}{A direction of the alternative hypothesis. A two-sided testing is assumed as a default (\code{two.sided}). Other choices include \code{greater} or \code{less}.}

\item{seed}{A random seed to be used for permutation. It is important to use the same seed for integrating results from both hemispheres. If \code{seed = NULL}, then a random number will be chosen.}

\item{parallel}{Whether parallel computing is to be used. Default: \code{FALSE}}

\item{cores}{The number of cores when parallel computing is executed.}
}

\description{
TBD
}
\examples{
TBD
}