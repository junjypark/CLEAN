# ymat        : a V times N matrix (V: # of vertices, N: # of images)
# mod0        : a N times p covariate matrix (p: # of covariates)
#             : It can be generated easily by using the model.matrix() function.
#             : Make sure the covariate of interest that will be tested
#             : is NOT included in mod0
# group       : 
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

Clean=function(ymat, 
               mod0 = NULL,
               group = NULL, 
               distmat = NULL, 
               sacf = "exp",
               max.radius = 20,
               nperm = 5000, 
               alpha = 0.05, 
               alternative = c("two.sided", "less", "greater"), 
               seed = NULL, 
               partition = T, 
               npartition = NULL, 
               parallel = F, 
               ncores = 1){
  
  if ( alpha < 0 |alpha > 1){
    stop("[CLEAN] alpha should range between 0 and 1.")
  }
  if ( nrows(distmat) != nrow(ymat) ){
    stop("[CLEAN] The number of rows of distmat and the number of rows of ymat needs to be the same (# vertices).")
  }
  if (length(alternative) > 1){ 
    alternative = "two.sided"
    cat("[CLEAN] Conducting the two-sided test as alternative has not been specified.\n")
  }
  if (is.null(seed)){ 
    seed = sample(1e6,1) 
  }
  
  if (is.null(group)){
    cat("[CLEAN] Testing for the mean as group is not specified. \n")
    cat("[CLEAN] mod0 will not be used. \n")
    
    fit = CleanMean(ymat = ymat, 
                    distmat = distmat, 
                    sacf = sacf,
                    nperm = nperm, 
                    alpha = alpha, 
                    alternative = alternative,
                    seed = seed,
                    partition = partition, 
                    npartition = npartition,
                    parallel = parallel,
                    ncores = ncores)
    set.seed(NULL)
    return(fit)
  }
  else{
    if (length(group) != ncol(ymat)){
      stop("[CLEAN] The number of elements in group does not match with the number of columns in ymat.")
    }
    
    fit = CleanDiff(ymat = ymat, 
                    mod0 = mod0,
                    group = group, 
                    distmat = distmat, 
                    sacf = sacf,
                    nperm = nperm, 
                    alpha = alpha, 
                    alternative = alternative,
                    seed = seed,
                    partition = partition, 
                    npartition = npartition,
                    parallel = parallel, 
                    ncores = ncores)
    set.seed(NULL)
    return(fit)
  }
}

