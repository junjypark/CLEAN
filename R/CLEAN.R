# ymat        : A V times N matrix (V: # of vertices, N: # of images)
# distmat     : A V times V distance matrix
# cortex      : A vector of vertex indices that will be included. 
#             : It is needed if your ymat and distmat contained medial walls that should not be included. 
#             : Default: NULL (all vertices will be used).
# mod0        : A N times p covariate matrix (p: # of covariates)
#             : It can be generated easily by using the model.matrix() function.
#             : Make sure the covariate of interest (cov.interest) that will be tested is NOT included in mod0
# cov.interest: A length N vector for the covariate of interest.
#             : It can be a binary character vector for two-sample testing.
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
               distmat = NULL, 
               cortex = NULL,
               mod0 = NULL,
               cov.interest = NULL, 
               sacf = "exp",
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
    seed = sample(1e6, 1) 
  }
  
  if (is.null(cov.interest)){
    cat("[CLEAN] Testing for the mean as cov.interest is not specified. \n")
    cat("[CLEAN] mod0 will not be used. \n")
    
    fit = CleanMean(ymat = ymat, 
                    distmat = distmat, 
                    cortex = cortex,
                    sacf = sacf,
                    nperm = nperm, 
                    alpha = alpha, 
                    alternative = alternative,
                    seed = seed,
                    nngp = nngp,
                    nngp.J = nngp.J,
                    partition = partition, 
                    npartition = npartition,
                    parallel = parallel,
                    ncores = ncores)
    set.seed(NULL)
    return(fit)
  }
  else{
    if (length(table(cov.interest)) == 2){
      cov.interest = ifelse(cov.interest == cov.interest[1], 1, -1)
    }
      
    if (length(cov.interest) != ncol(ymat)){
      stop("[CLEAN] The number of elements in cov.interest does not match with the number of columns in ymat.")
    }
    
    fit = CleanDiff(ymat = ymat, 
                    cortex = cortex,
                    mod0 = mod0,
                    cov.interest = cov.interest, 
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

