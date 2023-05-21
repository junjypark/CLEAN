# ymat        : a V times N matrix for the first modality (V: # of vertices, N: # of images)
# xmat        : a V times N matrix for the second modality (V: # of vertices, N: # of images)
# mod         : a N times p covariate matrix (p: # of covariates)
#             : It can be generated easily by using the model.matrix() function.
# distmat     : a V times V distance matrix
# sacf        : spatial autocorrelation function
#             : The exponential function is assumed as a default.
#             : Other choices include "gau" (Gaussian) 
#             : and "mix" (mixture of exponential and Gaussian)
# max.radius  : The maximum radius for cluster enhancement. 20 is assumed as a default.
# nperm       : number of permutations to be used. At least 5000 permutation is recommended.
# alpha       : A desired FWER. alpha=0.05 is assumed as a default.
# alternative : A direction of the alternative hypothesis. A two-sided testing is assumed as a default.
# seed        : A random seed to be used for permutation.
#             : It is important to use the same seed for integrating results from both hemispheres.
# parallel    : Whether parallel computing is to be used.
# cores       : The number of cores when parallel computing is executed.

CleanR=function(xmat, 
                ymat,
                distmat, 
                cortex = NULL,
                mod = NULL,
                sacf = "exp",
                max.radius = 20,
                nperm = 5000, 
                alpha = 0.05,
                alternative = c("two.sided", "less", "greater"),
                seed = NULL, 
                nngp = T,
                nngp.J = 50,
                npartition = NULL, 
                parallel = F, 
                ncores = 1){

  if ( alpha<0 | alpha>1){
    stop("alpha should range between 0 and 1.")
  }
  if (length(alternative) > 1){
    alternative = "two.sided"
    cat("Conducting the two-sided test as alternative has not been specified...\n")
  }
  if (is.null(seed)){ 
    seed=sample(1e6,1) 
  }


  V = nrow(ymat)
  if (!is.null(cortex)){
    ymat = ymat[cortex, ]
    xmat = xmat[cortex, ]
    distmat = distmat[cortex, cortex]
  } else{
    cortex= 1:V
  }
  
  if (is.null(mod)){
    mod=rep(1, ncol(ymat))
  }
  
  ymat.leverage = spLeverage(data=ymat, distmat=distmat, mod0=mod, sacf=sacf, nngp=nngp, nngp.J=nngp.J)$out
  xmat.leverage = spLeverage(data=xmat, distmat=distmat, mod0=mod, sacf=sacf, nngp=nngp, nngp.J=nngp.J)$out
  xmat.leverage=t(apply(xmat.leverage,1,scale))/sqrt(ncol(xmat.leverage)-1)
  ymat.leverage=t(apply(ymat.leverage,1,scale))/sqrt(ncol(ymat.leverage)-1)
  
  NNmatrix = buildNNmatrixDist(distmat, max.radius = max.radius)
  
  if (is.null(npartition)){
    npartition = nrow(NNmatrix)%/%10000+1
  }
  
  NNList = list()
  len=ceiling(nrow(NNmatrix)/npartition)
  for (i in 1:npartition){
    start=len*(i-1)+1
    end=min(len*i,nrow(NNmatrix))
    NNList[[i]] = NNmatrix[start:end,]
  }
  
  CleanerPerm = CleanerPermC(xmat.leverage, ymat.leverage, nperm, seed)
  if (isTRUE(parallel)){
    cl = makeCluster(ncores)
    registerDoParallel(cl)
    result = foreach(i = 1:npartition, .packages=("CLEAN"))%dopar%{
      fit = CleanerExpandPermC(CleanerPerm$U, CleanerPerm$permU, NNList[[i]])
      fit$alternative = alternative
      fit$seed = seed
      fit
    }
    stopCluster(cl)
  } else{
    result=list()
    for (i in 1:npartition){
      result[[i]] = CleanerExpandPermC(CleanerPerm$U, CleanerPerm$permU, NNList[[i]])
      result[[i]]$alternative = alternative
      result[[i]]$seed = seed
    }  
  }
  cat("temp1 \n")
  
  result = combine(result, alpha = alpha, collapse = T)
  cat("temp1 \n")
  result_proc = process(result)
  result$Tstat = rep(0, V)
  result$Tstat[cortex] = result_proc$Tstat
  result$Tstat_thresholded = rep(0, V)
  result$Tstat_thresholded[cortex] = result_proc$Tstat_thresholded

  set.seed(NULL)
  
  return(result)
}
