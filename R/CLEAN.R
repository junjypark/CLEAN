Clean=function(ymat, NNmatrix=NULL, group=NULL, nperm=10000, alpha=0.05, alternative=c("two.sided","less", "greater"), seed=NULL, 
               partition=T, npartition=NULL, parallel=F, ncores=1){
  if (is.null(NNmatrix)){
    cat("[CLEAN] As NNmatrix is not specified, Clean conducts massive univariate analysis.\n")
    if (is.null(group)){
      fit=MassiveMean(ymat=ymat,
                      nperm=nperm, 
                      alpha=alpha, 
                      alternative=alternative,
                      seed=seed,
                      partition=partition, 
                      npartition = npartition,
                      parallel=parallel,
                      ncores=ncores)
      set.seed(NULL)
      return(fit)
    } else{
      fit=MassiveDiff(ymat=ymat,
                      group=group,
                      nperm=nperm, 
                      alpha=alpha, 
                      alternative=alternative,
                      seed=seed,
                      partition=partition, 
                      npartition = npartition,
                      parallel=parallel,
                      ncores=ncores)
      set.seed(NULL)
      return(fit) 
    }
  } else{
    if (is.null(group)){
      fit=CleanMean(ymat=ymat, 
                    NNmatrix=NNmatrix, 
                    nperm=nperm, 
                    alpha=alpha, 
                    alternative=alternative,
                    seed=seed,
                    partition=partition, 
                    npartition = npartition,
                    parallel=parallel,
                    ncores=ncores)
      set.seed(NULL)
      return(fit)
    }
    else{
      fit=CleanDiff(ymat=ymat, 
                    NNmatrix=NNmatrix, 
                    group=group, 
                    nperm=nperm, 
                    alpha=alpha, 
                    alternative=alternative,
                    seed=seed,
                    partition=partition, 
                    npartition = npartition,
                    parallel=parallel, 
                    ncores=ncores)
      set.seed(NULL)
      return(fit)
    }
  }
}

