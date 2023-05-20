\name{combine}
\alias{combine}
\title{Combining results from different hemispheres/regions}
\usage{
combine(
  lst, 
  alpha=0.05,
  collapse = F
)
}

\arguments{

\item{lst}{The list of CLEAN (or CLEAN-R or CLEAN-V) objects across different hemispheres/regions. For each fit, \code{seed} and \code{nperm} must be the same.}

\item{alpha}{The 'brain-wise' FWER. Default: 0.05}

\item{collapse}{Whether to collapse the list into a single fit. Default: \code{F}}

}


\description{
Combining different CLEAN (or CLEAN-R or CLEAN-V) fits across different regions or hemispheres. This step is necessary to control FWER at the brain level. Therefore, after applying \code{combine()}, the \code{Tstat_thresholded} will become slightly more conservative and the \code{threshold} will be increased. 
}
