\name{Clean}
\alias{Clean}
\title{Cluster enhancement for testing a regression parameter in GLM}
\usage{
Clean(
  ymat,
  distmat = NULL,
  cortex = NULL,
  mod0 = NULL,
  group = NULL,
  sacf = "exp"
  max.radius = 20,
  nperm = 5000,
  alpha = 0.05,
  alternative = c("two.sided", "less", "greater"),
  seed = NULL,
  nngp = T,
  nngp.J = 50,
  partition = T,
  npartition = NULL,
  parallel = F,
  ncores = 1
)
}

\arguments{
\item{ymat}{A V times N matrix (V: # of vertices, N: # of images.}

\item{distmat}{A V times V distance matrix}

\item{cortex}{A vector of vertex indices that will be included. It is needed if your ymat and distmat contained medial walls that should not be included. Default: \code{NULL} (all vertices will be used).}

\item{mod0}{A N times p covariate matrix (p: # of covariates). It can be generated easily by using the \code{model.matrix()} function. Make sure the covariate of interest that will be testedis NOT included in mod0.}

\item{cov.interest}{A covariate vector of interest. For two-sample testing, it needs to be 1 for group 1 and -1 for group 2.}

\item{sacf}{A spatial autocorrelation function. The exponential function (\code{exp}) is assumed as a default. Other choices include \code{gau} (Gaussian) or \code{mix} (mixture of exponential and Gaussian).}

\item{max.radius}{The maximum radius for cluster enhancement. Default: 20.}

\item{nperm}{The number of permutations to be used. At least 5000 permutation is recommended. Default: 5000.}

\item{alpha}{A desired FWER. Default: 0.05.}

\item{alternative}{A direction of the alternative hypothesis. A two-sided testing is assumed as a default (\code{two.sided}). Other choices include \code{greater} or \code{less}.}

\item{seed}{A random seed to be used for permutation. It is important to use the same seed for integrating results from both hemispheres. If \code{seed = NULL}, then a random number will be chosen.}

\item{nngp}{Whether or not NNGP will be used for computing test statistics. It is useful for high dimensional data. Default: \code{TRUE}.}

\item{nngp.J}{The number of nearest neighbors used to operate NNGP (if \code{nngp = T}. It is different from \code{max.radius} and it is recommended not to change this value unless you are confident of its usage. Default: 50.}

\item{parallel}{Whether parallel computing is to be used. Default: \code{FALSE}}

\item{cores}{The number of cores when parallel computing is executed. Default: 1.}

}

\value{
A list of the following is returned.

\item{Tstat}{A length V vector of cluster enhanced test statistics obtained by CLEAN. Non-cortex vertices take the value 0.}

\item{Tstat_thresholded}{A cluster enhanced test statistics obtained by CLEAN, thresholded by a pre-set FWER level (alpha). Any vertices whose cluster-enhanced test statistics do not pass the threshold will be shurunk to 0. Non-cortex vertices take the value 0 as well. Use this with cautions because it would need to be adjusted for combining results from different hemispheres/regions. See \code{combine()} for details.}
}

\item{permMax}{A nperm vector that contains the maximum of cluster-enhanced test statistics across all vertices at each permuted sample. Used to compute the FWER-controlling threshold.}

\item{permMin}{A nperm vector that contains the minimum of cluster-enhanced test statistics across all vertices at each permuted sample. Used to compute the FWER-controlling threshold.}

\item{threshold}{A FWER-controlling threshold. Use this with cautions because it would need to be adjusted for combining results from different hemispheres/regions. See \code{combine()} for details.}

\item{nperm}{The number of permutation used.}

\item{seed}{The seed number used.}

\item{nlocations}{The number of cortical vertices used.}

\item{alternative}{The direction of the alternative hypothesis (\code{two.sided}, \code{greater}, or \code{less})}

\description{
Fitting CLEAN to test activation or brain-behavior associations in GLM.
}

\examples{
Forthcoming
}
