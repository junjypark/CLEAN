# Temporary version --- use with cautions

# ymat        : a V times N matrix of imaging data (V: # of vertices, N: # of images)
# mod         : a N times p covariate matrix (p: # of covariates)
#             : It can be generated easily by using the model.matrix() function.
# K           : a N times N matrix specifying between-image dependencies.
# distmat     : a V times V distance matrix
# sacf        : spatial autocorrelation function
#             : The exponential function is assumed as a default.
#             : Other choices include "gau" (Gaussian) 
#             : and "mix" (mixture of exponential and Gaussian)
# max.radius  : The maximum radius for cluster enhancement. 20 is assumed as a default.
# nperm       : number of permutations to be used. At least 5000 permutation is recommended.
# alpha       : A desired FWER. alpha=0.05 is assumed as a default.
# seed        : A random seed to be used for permutation.
#             : It is important to use the same seed for integrating results from both hemispheres.
# parallel    : Whether parallel computing is to be used.
# cores       : The number of cores when parallel computing is executed.

CleanV=function(ymat,
                distmat, 
                cortex = NULL,
                mod = NULL,
                K,
                sacf = "exp",
                max.radius = 20,
                nperm = 5000, 
                alpha = 0.05,
                seed = NULL, 
                nngp = T,
                nngp.J = 50,
                npartition = NULL, 
                parallel = F, 
                ncores = 1){
  
  if ( alpha < 0 |alpha > 1){
    stop("[CLEAN] alpha should range between 0 and 1.")
  }
  if ( nrows(distmat) != nrow(ymat) ){
    stop("[CLEAN] The number of rows of distmat and the number of rows of ymat needs to be the same (# vertices).")
  }
  if (is.null(seed)){ 
    seed = sample(1e6, 1) 
  }
  
  ymat.leverage = spLeverage(ymat, distmat, mod, sacf, nngp, nngp.J)
  NNmatrix = buildNNmatrixDist(distmat, max.radius = max.radius)
  
  K=Matrix(K, sparse=T)
  
  out = CleanVarC(ymat.leverage, NNmatrix, K, nperm, seed)
  
  out$seed = seed
  out$nlocations = ncol(NNmatrix)
  out$alternative = "greater"
  out = combine(out, alpha = alpha)
  
  set.seed(NULL)
  
  return(out)
}

